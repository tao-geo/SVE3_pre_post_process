import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# import xarray as xr

# for step in range(0, 81):
#     lon, lat, ocn = np.loadtxt(f'ICE6G_1x1_refine_ocean_iter2/ocn.{step}', skiprows=1, unpack=True)
#     fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
#     ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
#     cb = ax.scatter(lon, lat, c=ocn, s=0.1, cmap='viridis', transform=ccrs.PlateCarree())
#     ax.coastlines()
#     plt.colorbar(cb, ax=ax, orientation='horizontal', label='OCN function')
#     plt.title(f"Step {step}")
#     plt.savefig(f"./ICE6G_1x1_refine_ocean_iter2/ocn_{step}.png", dpi=200)
#     plt.close(fig)

for step in range(0, 81):
    lon, lat, ice = np.loadtxt(f'ICE6G_1x1_refine_ocean_iter2/ice.{step}', skiprows=1, unpack=True)
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    cb = ax.scatter(lon, lat, c=ice, s=0.1, cmap='viridis', transform=ccrs.PlateCarree())
    ax.coastlines()
    plt.colorbar(cb, ax=ax, orientation='horizontal', label='load')
    plt.title(f"Step {step}")
    plt.savefig(f"./ICE6G_1x1_refine_ocean_iter2/ice_{step}.png", dpi=200)
    plt.close(fig)