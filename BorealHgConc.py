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
from matplotlib.path import Path

# 设置全局字体为Arial
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

# HgConc = gdal.Open("BorealHg results/2021CC.tif").ReadAsArray()
# HgConc = gdal.Open("BorealHg results/2023CC.tif").ReadAsArray()
HgConc = gdal.Open("BorealHg results/15-19CC.tif").ReadAsArray()

Fillgap=HgConc[:,0]/2+HgConc[:,-1]/2
HgConc = np.column_stack((HgConc, Fillgap))

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection=ccrs.LambertAzimuthalEqualArea(central_latitude=90))
ax.set_extent([-180, 180, 45, 90], crs=ccrs.PlateCarree())

#裁剪圆形
theta = np.linspace(0, 2*np.pi, 100)
# 圆心 (0.5, 0.5)，半径 0.5，恰好内切于方形轴
verts = np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + 0.5
circle = Path(verts, closed=True)
ax.set_boundary(circle, transform=ax.transAxes)  # 应用圆形裁剪

lon = np.arange(-180, 180+2, 2.5)
lat = np.arange(90, -88-2, -2)
lon_grid, lat_grid = np.meshgrid(lon, lat)

bounds = np.arange(0, 0.275, 0.025)
# 颜色值（从左到右对应 -12 到 +12）
original_cmap = plt.get_cmap('YlOrRd')
colors = original_cmap(np.linspace(0, 0.8, 128))
cmap = ListedColormap(colors)

# 使用 contourf 并指定 levels
mesh = ax.contourf(lon_grid, lat_grid, HgConc, 
                   levels=bounds, 
                   cmap=cmap, 
                   transform=ccrs.PlateCarree(),
                   vmin=0, vmax=0.25,
                   extend='both')

ax.add_feature(cfeature.LAND, facecolor='none', edgecolor='black', linewidth=1, zorder=2)

# 北极圈
lon_circle = np.linspace(-180, 180, 200)
lat_arctic = np.full_like(lon_circle, 66.5)
ax.plot(lon_circle, lat_arctic,
        transform=ccrs.PlateCarree(),
        color='gray', dashes=(6, 4), linewidth=1, alpha=0.7)

# --- 手动添加经度标签 ---
def format_lon(lon):
    if lon == 0:
        return "0°"
    elif lon == 180:
        return "180°"
    elif lon > 0:
        return f"{int(lon)}°E"
    else:
        return f"{int(-lon)}°W"

target_lons = np.arange(-150, 180+30, 30)   # 每隔30°一条经线
label_lat_dict = {
    0: 42,
    180: 42,
    60: 38,
    -60: 38,
    120: 38,
    -120: 38,
}

for lon0 in target_lons:
    if lon0 % 60 == 0:
        label_lat = label_lat_dict.get(lon0, 40)
        ax.text(lon0, label_lat, format_lon(lon0),
                transform=ccrs.PlateCarree(),
                ha='center', va='center',
                fontsize=14, color='black')

# 颜色条
cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', 
                    pad=0.06, shrink=0.5, aspect=30, extend='both',
                    ticks=np.arange(0, 0.255, 0.05))
cbar.ax.tick_params(labelsize=10)

plt.title('Hg$^{0}$ concentration (ng m$^{-3}$)',fontsize=12, pad=20, y=-0.26, fontname='Arial')
# plt.savefig('BorealHg results/2021CCfromMBOBB.png', transparent=True, dpi=300)
# plt.savefig('BorealHg results/2023CCfromMBOBB.png', transparent=True, dpi=300)
plt.savefig('BorealHg results/15-19CCfromMBOBB.png', transparent=True, dpi=300)
plt.show()
