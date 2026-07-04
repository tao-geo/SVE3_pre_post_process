
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


def plot_on_map(lon, lat, data, fig_text='analytic', cmap='coolwarm', vcenter=None, projection='P', save_filename=None):
    '''
    parameters:
    1. lon, lat, data should be (n,1) array
    2. if vcenter is set, the colormap is a two slope norm
    3. projection: 'P' or 'M'
    4. save_filename: if set, figure is saved as the argument passed
    '''

    fig = plt.figure(figsize=(16, 10))
    if projection=='M':
        proj = ccrs.Mollweide()
    elif projection=='P':
        proj = ccrs.PlateCarree(center_longitude=180)
    elif projection=='Robinson':
        proj = ccrs.Robinson(central_longitude=180)
    elif projection=='AlbersEqualArea':
        proj = ccrs.AlbersEqualArea()
    elif projection=='LambertConformal':
        proj = ccrs.LambertConformal(cutoff=0)
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.coastlines()
    # Add grid lines and tick labels to the map
    ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
    cmap = 'coolwarm'
    if vcenter!=None:
        from matplotlib.colors import TwoSlopeNorm
        norm = TwoSlopeNorm(vmin=np.min(data), vcenter=0, vmax=np.max(data))
        cb0 = ax.scatter(lon, lat, c = data, transform=ccrs.PlateCarree(),cmap=cmap, norm=norm) #, cmap='Greys' )
        # cb0 = ax.contourf(lon,lat,data,transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
    else:
        cb0 = ax.scatter(lon, lat, c = data,s=1 , transform=ccrs.PlateCarree(),cmap=cmap ) 
    plt.colorbar(cb0, ax=ax, location='bottom', fraction=0.1, aspect=30)
    ax.title.set_text(fig_text)
    if save_filename != None:
        plt.savefig(save_filename)

def remove_duplicate(lon, lat, visc):
    '''
    Remove duplicate points in the data.
    '''
    lonlat = np.vstack([lon, lat]).T
    _, idx = np.unique(lonlat, axis=0, return_index=True)
    lon = lon[idx]
    lat = lat[idx]
    visc = visc[idx]
    return (lon, lat, visc)

#################### main ####################

## I. read data

depth_id = 3 #4
depth= 30 #60
case_id = 'case8_WPM21'

data = np.loadtxt(f'../data_for_download/{case_id}/{case_id}.map_visc.{depth_id}.1')
lon = data[:,0]
lat = data[:,1]
visc = data[:,2]

## II. remove duplicate points
lon, lat, visc = remove_duplicate(lon, lat, visc)
visc = np.log10(visc)
## III. plot
fn_fig = f'../data_for_download/{case_id}/visc_{depth}km.png'
plot_on_map(lon, lat, visc, fig_text=f'viscosity at {depth}km', cmap='coolwarm', 
             projection='Robinson', save_filename=fn_fig)
print(f'Figure saved as {fn_fig}')
