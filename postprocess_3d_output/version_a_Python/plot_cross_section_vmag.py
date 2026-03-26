'''
plot V_mag (Vr^2 + Vt^2 + Vphi^2)^1/2 along certain cross section.
'''


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pygmt



steps = np.arange(16, 321, 16)

R = 6371e3



steps = np.arange(16, 321, 16)

R = 6371e3

lon2plot = 250
# # lon2plot = 
# # great circle along a path
# p1 = [220, 40]
# p2 = [320, 40]

case_id = 'case_A_ICE0_V1D'
WITH_ASTHENO = False

# case_id = 'case_B_ICE0'
# WITH_ASTHENO = True

fn_fig_out = f'../figs/{case_id}_Vmag_slice_lon{lon2plot}'


# ice 
fn_ice = '/glade/work/taoyuan/Proj_RSL_3Dvisc/CASES/preprocess_initial_40ka_NEW/IceT.I6F_C.131QB_VM5a_1deg.nc'
xr_ice = xr.open_dataset(fn_ice)
time_ice = xr_ice['Time'].values # positive values
xr_ice_ih = xr_ice['stgit'] # time, lat, lon
# xr_ice_ref = xr_ice_ih[np.argmin(np.abs(time_ice-40)),:,:]

for i, step in enumerate(steps):
    fig = pygmt.Figure()

    # data related to time
    fn_grid = f'/glade/work/taoyuan/Proj_RSL_3Dvisc/CASES/postp_3D_disp/{case_id}_mantleOutput/{case_id}_mantleOutput.{step}.grid_lon200-330_lat0-90_res1.0.nc'
    ds = xr.open_dataset(fn_grid)
    years = ds.attrs['time'] / 1e3 - 40
    dts = ds.attrs['dt']

    # plot ice 
    fig.basemap(region=[220, 320, 0, 90], projection="Q6c", frame=True)
    pygmt.makecpt(cmap='vik', series="-100/100/2", reverse=True)
    itime_ice = np.argmin(np.abs(time_ice-np.abs(years)))
    ice_rate = (xr_ice_ih[itime_ice,:,:] - xr_ice_ih[itime_ice-2,:,:])/np.abs(time_ice[itime_ice]-time_ice[itime_ice-2])
    fig.grdimage(ice_rate)
    fig.coast(shorelines='0.1p,black', area_thresh=2e4)
    fig.colorbar(frame='a100f+l"ice thickness change rate (m/ka)"')
    # add section line
    fig.plot(x=[lon2plot, lon2plot], y=[0, 90], pen='1p,black')

    fig.shift_origin(xshift="7c")

    ######## plot V



    lat = ds['lat'][:]
    lon = ds['lon'][:]
    r = ds['r'][:]

    ilon_selected = np.argmin(np.abs(lon.values-lon2plot))
    vr = ds['dr_incr'][:,:,ilon_selected]*R/dts * 100 # cm/yr
    vtheta = ds['dtheta_incr'][:,:,ilon_selected]*R/dts * 100 # cm/yr
    vphi = ds['dphi_incr'][:,:,ilon_selected]*R/dts * 100 # cm/yr


    dr_selected = np.sqrt(vr**2 + vtheta**2 + vphi**2)

    dr_selected = dr_selected.drop_duplicates(dim='r')

    r_equal_space = np.linspace(np.min(r.values), np.max(r.values), 100)
    dr_equal_space = dr_selected.interp(r=r_equal_space)

    fig.basemap(region=[0, 90, np.min(dr_equal_space['r'].values), np.max(dr_equal_space['r'].values)], projection='X12c/6c', 
                frame=[f'WSne+tt={years:.1f} BP', 'xaf+lLat', 'yaf+lr'])

    vmax = None

    if np.abs(years - -14) < 0.1:
        vmax = 10
    elif years < -29:
        vmax = 0.5
    elif years < -6:
        vmax = 5
    else:
        vmax = 2

    pygmt.makecpt(cmap='magma', series=[0, 2*vmax], reverse=True)
    fig.grdimage(dr_equal_space)
    # fig.contour(dr_equal_space.values, levels=1, pen='0.5p')

    # plot Asthenosphere and 670,
    if WITH_ASTHENO:
        fig.plot(x=[0, 90], y=[1-300e3/R, 1-300e3/R], pen='1p,black')
        fig.text(x=0, y=1-300e3/R, text='Asthenosphere', font='6p,Helvetica-Bold,black', justify='LB', offset='5p/2p')
    else:
        fig.plot(x=[0, 90], y=[1-300e3/R, 1-300e3/R], pen='1p,black,--')
        fig.text(x=0, y=1-300e3/R, text='300km', font='6p,Helvetica-Bold,black', justify='LB', offset='5p/2p')
    
    fig.plot(x=[0, 90], y=[1-670e3/R, 1-670e3/R], pen='1p,red')
    fig.text(x=0, y=1-670e3/R, text='670 km', font='6p,Helvetica-Bold,red', justify='LB', offset='5p/2p')
    fig.text(x=0, y=np.min(dr_equal_space['r'].values), text='CMB', font='6p,Helvetica-Bold,red', justify='LB', offset='5p/2p')
    fig.colorbar(frame='af+l"velo magnitude(cm/yr)"')
    # if i%10 == 0:
    # fig.show()
    fig.savefig(f'{fn_fig_out}.{i}.png')



###########################


