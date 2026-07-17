import numpy as np
import pygmt
import shtns
# import gravity_toolkit




def gen_xarray(lat, lon, data):
    import xarray as xr
    xr_data = xr.DataArray(
        data,
        coords=[lat, lon],
        dims=["lat", "lon"],
    )
    return xr_data

def sum_sph(coeff, lmax):
    mmax = lmax
    sh = shtns.sht(lmax, mmax)  # create sht object with given lmax and mmax (orthonormalized)
    nlat, nphi =  (4*(lmax+1), 4*2*(lmax+1))
    sh.set_grid(nlat=nlat, nphi=nphi)  # build grid (gauss grid, phi-contiguous

    lat_grid_out = np.arcsin(sh.cos_theta) * 180 / np.pi
    # lat_grid_out = lat_grid_out[::-1] # from south pole to north pole
        # Note: SHTns require lat from north pole to south pole, here we reverse it.
        #   so, later, when we calculate the spherical harmonics, we need to reverse the latitudes.

    lon_grid_out = (2.*np.pi/nphi)*np.arange(nphi) * 180 / np.pi

    # data format is (l, m, clm, slm), i.e. degree, order, geoid stokes coeff, uplift Stokes, and horizontal Stokes

    ylm = sh.spec_array()
    gylm = sh.spec_array()

    G = 6.67430e-11  # m^3 kg^-1 s^-2
    factor0 = 4 * np.pi * G * 4400 * (6371e3)**2
    i=0
    for l in range(lmax+1):
        for m in range(l+1):
            if m == 0:
                factor = factor0
            else:
                factor = factor0 / np.sqrt(2)
            if (l!=coeff[i,0]) or (m!=coeff[i,1]):
                print("something is wrong!", l, m, coeff[i,0], coeff[i,1])
                exit()
            idx = sh.idx(l, m)
            ylm[idx] = coeff[i,2]* factor + 1j*coeff[i,3] * (-factor)
            gylm[idx] = (l - 1)/6371e3 * ylm[idx]

            i = i+1
    ylm /= 9.82
    geoid = sh.synth(ylm)
    gravity_anomaly = sh.synth(gylm)

    xr_geoid = gen_xarray(lat_grid_out, lon_grid_out, geoid)
    xr_gravity_anomaly = gen_xarray(lat_grid_out, lon_grid_out, gravity_anomaly)

    return xr_geoid, xr_gravity_anomaly

################################################# calculation

case_name = 'case_A_ICE0_V1D_mantleOutput'
case_dir = '/glade/work/taoyuan/Proj_RSL_3Dvisc/CASES/postprocess/'
fn_gdot = f"{case_dir}/copied_{case_name}/{case_name}.pttldot_sharm.320"
dt = 125 # yr
coeff = np.loadtxt(fn_gdot, skiprows=2)
coeff[:,2:] /= dt

fn_geoid_rate_from_file = '/glade/work/taoyuan/Proj_RSL_3Dvisc/CASES/postprocess/combined_files_case_A_ICE0_V1D/case_A_ICE0_V1D.map_rate_geoid.320.regular.grd'

xr_geoid, xr_gravity_anomaly = sum_sph(coeff, 32)


############################ check 

fig = pygmt.Figure()

fig.basemap(region="g", projection='Q6i', frame=['WSne+tGeoid rate from sph'])
fig.grdimage(xr_geoid, cmap='vik')
fig.colorbar()
fig.coast(shorelines=True, area_thresh=2e5)

fig.shift_origin(xshift="6.5i")
fig.basemap(region="g", projection='Q6i', frame=['WSne+tGeoid rate from file'])
fig.grdimage(fn_geoid_rate_from_file, cmap='vik')
fig.colorbar()
fig.coast(shorelines=True, area_thresh=2e5)

import xarray as xr
data_from_file = xr.load_dataset(fn_geoid_rate_from_file)['z']
xr_diff = data_from_file - xr_geoid.interp_like(data_from_file)
fig.shift_origin(xshift="6.5i")
fig.basemap(region="g", projection='Q6i', frame=['WSne+tGeoid rate diff'])
fig.grdimage(xr_diff, cmap='vik')
fig.colorbar()
fig.coast(shorelines=True, area_thresh=2e5)


fig.shift_origin(xshift="6.5i")
fig.basemap(region="g", projection='Q6i', frame=['WSne+tGravity Anomaly'])
fig.grdimage(xr_gravity_anomaly, cmap='vik')
fig.colorbar()
fig.coast(shorelines=True, area_thresh=2e5)

fig.savefig(f"test/geoid_{case_name}.png")