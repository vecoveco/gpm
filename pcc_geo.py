import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import wradlib as wrl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
'''
def plot_dem(ax):
    filename = wrl.util.get_wradlib_data_file('geo/Europe_2_02.2017075.aqua.ndvi.250m.tif')
    # pixel_spacing is in output units (lonlat)
    rastercoords, rastervalues = wrl.io.read_raster_data(filename,
                                                         spacing=0.005)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                        rastervalues, cmap=plt.cm.terrain, norm=LogNorm(),
                        vmin=1, vmax=3000)
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    cb = plt.gcf().colorbar(dem, cax=cax1,
                           ticks=ticker.LogLocator(subs=range(10)))
    cb.set_label('terrain height [m]')
'''
def plot_dem(ax):
    filename = wrl.util.get_wradlib_data_file('geo/radolan_900x900_cr_500.tif')
    # pixel_spacing is in output units (lonlat)
    rastercoords, rastervalues = wrl.io.read_raster_data(filename)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                        rastervalues+1, cmap=plt.cm.terrain, norm=LogNorm(),
                        vmin=1, vmax=3000)
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    cb = plt.gcf().colorbar(dem, cax=cax1,
                           ticks=ticker.LogLocator(subs=range(10)))
    cb.set_label('terrain height [m]')

from pcc import plot_borders
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, aspect='equal')
plot_dem(ax)
plot_borders(ax)

plt.show()

'''
filename1 = wrl.util.get_wradlib_data_file('geo/radolan_900x900_cr_500.tif')
# pixel_spacing is in output units (lonlat)
rastercoords, rastervalues = wrl.io.read_raster_data(filename1)#,
                                                     #spacing=0.005)



plt.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                       rastervalues)
plt.colorbar()
plt.show()
'''
#plt.scatter (h_terra, h_aqua)
#plt.grid()
#plt.xlabel('Terra Hights in m')
#plt.ylabel('Aqua Hights in m')
#plt.show()
