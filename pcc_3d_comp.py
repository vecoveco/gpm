"""

3d Radar composit

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr
from pcc import get_time_of_gpm
from pcc import cut_the_swath
from netCDF4 import Dataset
from wradlib.io import read_generic_netcdf
from wradlib.util import get_wradlib_data_file
import os
## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders
from pcc import plot_radar
from pcc import get_miub_cmap
my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()

def read_and_overview(filename):
    """Read NetCDF using read_generic_netcdf and print upper level dictionary keys
    """
    test = read_generic_netcdf(filename)
    print("\nPrint keys for file %s" % os.path.basename(filename))
    for key in test.keys():
        print("\t%s" % key)

pfad = ('/automount/ags/velibor/gpmdata/nicht3dComposit/*.nc')

pfad_3d = sorted(glob.glob(pfad))[1]

comp3d = get_wradlib_data_file(pfad_3d)
#read_and_overview(comp3d)

from netCDF4 import Dataset
import numpy as np

comp3d = Dataset(pfad_3d, mode='r')

x = comp3d.variables['x'][:][300:600]
y = comp3d.variables['y'][:][300:600]
z = comp3d.variables['z'][:][1:60]

zh = comp3d['cband_oase_zh'][1:60,300:600,300:600]

#rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5
'''
from pcc import plot_borders
for i in range(len(x)):
    print str(i) + ' von ' + str(len(x))
    fig = plt.figure(figsize=(12,8))
    ax2 = fig.add_subplot(121, aspect='equal')
    plt.pcolormesh(x,y,zh[3,:,:], vmin=0, vmax=50)
    #plt.plot(x[i],y[i])
    #plt.axvline(x=i)
    plt.axline(x=x[i], color='k', linestyle='--')
    plt.colorbar(shrink=0.5)
    plt.grid()
    plot_borders(ax2)

    ax23 = fig.add_subplot(122, aspect='auto')
    plt.pcolormesh(y,z,zh[:,i,:], vmin=0, vmax=50)
    plt.colorbar(shrink=0.5)
    plt.grid()
    #plot_borders(ax23)

    #plt.title('Hight' + str(x[i]) + ' km')
    #plt.show()
    plt.tight_layout()
    plt.savefig('/automount/ags/velibor/plot/nicht3Dcomp/comp3d_x_'+str(x[i]) + '.png' )
    plt.close()
'''
'''
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x, y, z = x[::10], y[::10], z[::6]
zh = zh[::6, ::10, ::10]

#X, Y, Z = np.meshgrid(x, y, z)# (900, 900, 60)
#Z, Y, X = np.meshgrid(z, y, x)# (900, 60, 900)
Y, Z, X = np.meshgrid(y, z, x)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.scatter(X,Y,Z, c=zh)
plt.show()
'''


import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

X,Y,Z = np.meshgrid(x, y, z)
cc = zh

cc[cc<50]=np.nan
cc[cc==76.20000458]=np.nan


fig = plt.figure(1)
fig.clf()
ax = Axes3D(fig)

comp = ax.scatter(X,Y,Z, c=cc, alpha=0.3, cmap=my_cmap)

fig.colorbar(comp, shrink=0.5, aspect=5)

plt.draw()
plt.show()