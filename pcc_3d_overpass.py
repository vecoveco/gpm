"""

Einlesen und darstellen von GPM DPR und Radolan Dateien

Radolanpfad:

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

import csv

# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 12#18#7

pfad = ('/automount/ags/velibor/gpmdata/dpr/*.HDF5')
pfad_gpm = sorted(glob.glob(pfad))

print pfad

## Read GPM Data
## -------------
#try:
pfad_gpm_g = pfad_gpm[34]

print pfad_gpm_g

gpmdpr = h5py.File(pfad_gpm_g, 'r')

gprof_lat = np.array(gpmdpr['NS']['Latitude'])
gprof_lon = np.array(gpmdpr['NS']['Longitude'])
gpm_z = np.arange(0,176,1)

gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrected'])
gprof_pp[gprof_pp==-9999.9]= np.nan

print gprof_lat.shape, gprof_lon.shape, gprof_pp.shape

## Cut the GPM Swath
## ------------------


blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp,eu=0)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

print gpm_x.shape, gpm_z.shape, gprof_pp_b[:,25,:].shape

for i in range(0,160,5):
    plt.subplot(1,2,1)
    plt.pcolormesh(gpm_z,gpm_x[:,25], np.ma.masked_invalid(gprof_pp_b[:,25,:]), vmin=0, vmax=50)
    plt.subplot(1,2,2)
    plt.pcolormesh(gpm_x,gpm_y, np.ma.masked_invalid(gprof_pp_b[:,:,i]), vmin=0, vmax=50)
plt.show()

gx = gpm_x.ravel()
gy = gpm_y.ravel()
gz = gpm_z.ravel()

X,Y,Z = np.meshgrid(gx, gy, gz)
'''
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

#gx = gpm_x.reshape(gpm_x.shape[0]*gpm_x.shape[1])
#gy = gpm_y.reshape(gpm_y.shape[0]*gpm_y.shape[1])
gz = gpm_z

gx = gpm_x.ravel()
gy = gpm_y.ravel()

X,Y,Z = np.meshgrid(gx, gy, gz)
cc = gprof_pp_b

print X.shape, Y.shape, Z.shape, cc.shape

cc[cc<20]=np.nan

fig = plt.figure(1)
fig.clf()
ax = Axes3D(fig)

comp = ax.scatter(X,Y,Z, c=cc, alpha=0.3, cmap=my_cmap)

fig.colorbar(comp, shrink=0.5, aspect=5)

plt.draw()
plt.show()
'''