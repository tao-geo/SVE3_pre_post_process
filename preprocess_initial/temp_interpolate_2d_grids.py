# temporarily use python to interpolate among 2d regular grids (maybe not equal space)

import numpy as np
from scipy.interpolate import RegularGridInterpolator

PRESENT_TOPO_FILE = './ETOPO1+bedmap2_256.xyz'
FN_SAVE_TOPO = './ETOPO1+bedmap2_180x360.txt'

def main():


    # get grid
    nlat = 180
    nphi = 360

    # lat_grid_out = np.arcsin(sh.cos_theta) * 180 / np.pi
    # lat_grid_out = lat_grid_out[::-1] # from south pole to north pole

    lat_grid_out = np.linspace(-89.5, 89.5, nlat)
    # print("lat grid: ", lat_grid_out)

    lon_grid_out = np.linspace(0.5, 359.5, nphi)

    # load_present_topo(lon_grid_out, lat_grid_out)

    plot_topo()

    return

def load_present_topo(lon_out, lat_out):
    ''' read present day topo; interpolate to out grid
    Input: 
        lon_out, lat_out: 1D array:
    Return:
        topo_out: 2D array, (lat, lon), south pole to north pole
    '''
    fn_topo = PRESENT_TOPO_FILE

    topo = np.loadtxt(fn_topo)[:,-1]
    grid = np.loadtxt(fn_topo)[:,:-1]

    lon = grid[:,0]; lat = grid[:,1]

    topo_pd = topo.reshape((256,256*2)).T  # lon, lat
    lon_1d = grid[:,0].reshape((256,256*2))[0,:]
    lat_1d = grid[:,1].reshape((256,256*2))[:,0]

    lat_1d = lat_1d[::-1] # from south pole to north pole
    topo_pd = topo_pd[:,::-1] 

    # project to our grid
    
    # topo_out_func = RegularGridInterpolator((lon_1d, lat_1d), topo_pd, method="nearest", bounds_error=False, fill_value=None)
    topo_interp_func = RegularGridInterpolator((lon_1d, lat_1d), topo_pd, method='nearest', bounds_error=False, fill_value=None)
    lon_out_2d, lat_out_2d = np.meshgrid(lon_out, lat_out, indexing='ij')
    topo_out = topo_interp_func((lon_out_2d, lat_out_2d))
    topo_out = topo_out.reshape((len(lon_out), len(lat_out)))

    # write to file (lon, lat, topo_1d)
    lat_out_temp = lat_out[::-1]
    topo_out_temp = topo_out[:,::-1] # revert latitude
    with open(FN_SAVE_TOPO, 'w') as fout:
        for i in range(len(lat_out_temp)):
            for j in range(len(lon_out)):
                lon_temp = lon_out[j]
                lat_temp = lat_out_temp[i]
                data = topo_out_temp[j,i]
                fout.write(f"{lon_temp:.5f} {lat_temp:.5f} {data:.5e}\n")

    return #topo_out.T


def plot_topo():
    import matplotlib.pyplot as plt
    topo = np.loadtxt(FN_SAVE_TOPO)
    lon = topo[:,0]
    lat = topo[:,1]
    data = topo[:,2]
    c = plt.scatter(lon, lat, c=data, s=0.1)
    plt.colorbar()
    plt.savefig(f"{FN_SAVE_TOPO}.png")



main()