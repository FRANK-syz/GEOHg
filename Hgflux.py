import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.cm import get_cmap
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER

# 设置全局字体为Arial
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

Hg_flux = gdal.Open("OBBHg results/MultiYear_Total_flux.tif").ReadAsArray()
Hg_flux = np.where(Hg_flux <= 0, np.NaN, Hg_flux)
Hg_flux = Hg_flux[0:600, :]

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.OCEAN, facecolor='#edfbff', alpha=0.8)
ax.add_feature(cfeature.LAND, facecolor='#f0f0f0')  # 陆地浅灰色


lon = np.arange(-180, 180, 0.25)
lat = np.arange(90, -60, -0.25)
lon_grid, lat_grid = np.meshgrid(lon, lat)

# 使用RdYlBu颜色映射并反转
cmap = get_cmap('RdYlBu_r')  # 添加'_r'后缀来反转颜色条

# 手动设置分类边界（示例）
bounds = [0, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
n_classes = len(bounds) - 1
norm = BoundaryNorm(bounds, cmap.N)

# 使用新的颜色映射和标准化
mesh = ax.pcolormesh(lon_grid, lat_grid, Hg_flux, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.3)  # 海岸线

# 添加网格线
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0, linestyle='--')
gl.top_labels = False    # 不显示顶部标签
gl.right_labels = False  # 不显示右侧标签
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.ylocator = mticker.FixedLocator(np.arange(-60, 91, 30))

cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', 
                    pad=0.08, shrink=0.9, aspect=40,extend="both",
                    ticks=bounds,  # 显示所有边界值
                    format=mticker.ScalarFormatter())  # 使用标量格式化器
cbar.ax.tick_params(labelsize=14)

# 设置颜色条标签
tick_labels = [str(b) for b in bounds]
cbar.ax.set_xticklabels(tick_labels)

plt.title(r'Emission intensity (μg m$^{-2}$ yr$^{-1}$)', fontsize=18, pad=20, y=-0.42, fontname='Arial')
plt.savefig('OBBHg results/MultiYear_Total_flux.png', transparent=True, dpi=300)
plt.show()
