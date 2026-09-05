import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

# 设置全局字体为Arial
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

Hg0Dep = gdal.Open("BorealHg results/GFED_Hg0dep.tif").ReadAsArray()
# Hg0Dep = gdal.Open("BorealHg results/MBOBB_Hg0dep.tif").ReadAsArray()

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-178.75, 178.75, -89, 89], crs=ccrs.PlateCarree())


lon = np.arange(-178.75, 178.75+2.5, 2.5)
lat = np.arange(89, -89-2, -2)
lon_grid, lat_grid = np.meshgrid(lon, lat)

bounds = np.arange(0, 5+0.5, 0.5)
original_cmap = plt.get_cmap('YlOrRd')
colors = original_cmap(np.linspace(0, 0.8, 32))
cmap = ListedColormap(colors)

# 使用 contourf 并指定 levels
mesh = ax.contourf(lon_grid, lat_grid, Hg0Dep, 
                   levels=bounds, 
                   cmap=cmap, 
                   transform=ccrs.PlateCarree(),
                   vmin=0, vmax=5,
                   extend='both')

ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black', linewidth=0.5, zorder=2)

# 地图网格
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 30))

# 颜色条
cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', 
                    pad=0.05, shrink=0.7, aspect=40, extend='both',
                    ticks=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
cbar.ax.tick_params(labelsize=10)

plt.title('Hg$^{0}$ deposition (μg m$^{-2}$ yr$^{-1}$)',fontsize=12, pad=20, y=-0.26, fontname='Arial')
plt.savefig('BorealHg results/GFEDHg0Dep.png', transparent=True, dpi=300)
# plt.savefig('BorealHg results/MBOBBHg0Dep.png', transparent=True, dpi=300)

plt.show()
