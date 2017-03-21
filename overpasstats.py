"""

Program fuer Statistik er Overpasses Daten

"""

import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

a = pd.read_csv('/automount/user/velibor/SHKGPM/prog/output3radolandpr.csv', sep=',')


fr = [1,2,3]
tit = ['Maximal detected Rainrate','Minimal detected Rainrate','Mean Rainrate']

for ii in fr:

    plt.subplot(3,1,ii)
    plt.plot(a.index, a[[ii]])
    plt.title(tit[ii-1])
    plt.ylabel('Rainrate in mm/h')
    plt.xlim(0, a.shape[0])
    plt.grid()


plt.xlabel('timesteps of overpasses')
plt.show()


import wradlib as wrl
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker

def plot_dem(ax):
    filename = wrl.util.get_wradlib_data_file('geo/bonn_new.tif')
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

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, aspect='equal')
plot_dem(ax)
plt.show()