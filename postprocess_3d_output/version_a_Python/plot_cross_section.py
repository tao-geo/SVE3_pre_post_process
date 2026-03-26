'''
plot Vr along certain cross section.

The cross section should be along meridian (constant latitude)
'''


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pygmt



steps = np.arange(16, 321, 16)

R = 6371e3

lon2plot = 270
# # lon2plot = 
# # great circle along a path
# p1 = [220, 40]
# p2 = [320, 40]

# case_id = 'case_A_ICE0_V1D'
# title='VM5a-1D'
# WITH_ASTHENO = False

# case_id = 'case_B_ICE0'
# title='MZ21-3D'
# WITH_ASTHENO = True

case_id = 'case_L17_ICE6G_V1D'
title='L17-1D'
WITH_ASTHENO = False

fn_fig_out = f'../figs/{case_id}_Vr_slice_lon{lon2plot}'


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

    ######## plot Vr



    lat = ds['lat'][:]
    lon = ds['lon'][:]
    r = ds['r'][:]


    ilon_selected = np.argmin(np.abs(lon.values-lon2plot))
    Vr_selected = ds['dr_incr'][:,:,ilon_selected]*R/dts * 100 # cm/yr

    Vr_selected = Vr_selected.drop_duplicates(dim='r')

    r_equal_space = np.linspace(np.min(r.values), np.max(r.values), 100)
    Vr_equal_space = Vr_selected.interp(r=r_equal_space)

    fig.basemap(region=[0, 90, np.min(Vr_equal_space['r'].values), np.max(Vr_equal_space['r'].values)], projection='X12c/6c', 
                frame=[f'WSne+t{years:.1f} ka ({title})', 'xaf+lLatitude', 'yaf+lr'])

    vmax = None

    if np.abs(years - -14) < 0.1:
        vmax = 10
    elif years < -29:
        vmax = 0.5
    elif years < -6:
        vmax = 5
    else:
        vmax = 2

    pygmt.makecpt(cmap='vik', series=[-vmax, vmax])
    fig.grdimage(Vr_equal_space)
    # fig.contour(dr_equal_space.values, levels=1, pen='0.5p')

    #### plot horizontal velocity
    V_x = -ds['dtheta_incr'][:,:,ilon_selected]*R/dts * 100 # cm/yr
    V_x = V_x.drop_duplicates(dim='r')
    V_x = V_x.interp(r=r_equal_space)

    V_y = Vr_equal_space * 0.0

    # scale for vectors (x / per cm)
    scale = 0.5
    if years > -15:
        scale = 2

    # use grdvector by calling C API
    with pygmt.clib.Session() as ses:
        with ses.virtualfile_from_grid(V_x) as vx, ses.virtualfile_from_grid(V_y) as vy:
            args = f"{vx} {vy} -Ix10/10 -Q0.2c+e -S{scale}c+s0.5 -l5mm/yr -W1p,blue"
            ses.call_module("grdvector", args)

    fig.legend(position='BL')

    # plot Asthenosphere and 670,
    if WITH_ASTHENO:
        fig.plot(x=[0, 90], y=[1-300e3/R, 1-300e3/R], pen='1p,black')
        fig.text(x=0, y=1-300e3/R, text='Asthenosphere', font='6p,Helvetica-Bold,black', justify='LB', offset='5p/2p')
    else:
        fig.plot(x=[0, 90], y=[1-300e3/R, 1-300e3/R], pen='1p,black,--')
        fig.text(x=0, y=1-300e3/R, text='300km', font='6p,Helvetica-Bold,black', justify='LB', offset='5p/2p')
    
    fig.plot(x=[0, 90], y=[1-670e3/R, 1-670e3/R], pen='1p,red')
    fig.text(x=0, y=1-670e3/R, text='670 km', font='6p,Helvetica-Bold,red', justify='LB', offset='5p/2p')
    fig.text(x=0, y=np.min(Vr_equal_space['r'].values), text='CMB', font='6p,Helvetica-Bold,red', justify='LB', offset='5p/2p')
    fig.colorbar(frame='af+l"Vr (cm/yr)"')
    # if i%10 == 0:
    # fig.show()
    fig.savefig(f'{fn_fig_out}.{i}.png')



###########################


