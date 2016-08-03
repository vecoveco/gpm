"""Einlesen und darstellen von GPM IMERG Dateien"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib as mpl

def plot_borders(ax):

    from osgeo import osr
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)

    filename = wradlib.util.get_wradlib_data_file('geo/ne_10m_admin_0_'
                                              'countries.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')
    wradlib.vis.add_lines(ax, borders, color='black', lw=1, zorder=4)
    ax.autoscale()

## Einlesen von IMERG GPM global
## -----------------------------

pfad = ('/user/velibor/SHKGPM/data/imerg/*.HDF5')
pfad_imerg = sorted(glob.glob(pfad))
pfad_imerg = pfad_imerg[0]

gpmi = h5py.File(pfad_imerg,'r')
print ('IMERGE Variablen:', gpmi[u'Grid'].keys())

gpmi_lat = gpmi[u'Grid'][u'lat']
gpmi_lon = gpmi[u'Grid'][u'lon']
#gpmi_pre = gpmi[u'Grid'][u'HQprecipitation']
#gpmi_pre = gpmi[u'Grid'][u'precipitationUncal']
#gpmi_pre = gpmi[u'Grid'][u'probabilityLiquidPrecipitation']
gpmi_pre = gpmi[u'Grid'][u'IRprecipitation']


gpmi_pre = np.array(gpmi_pre).transpose()
gpmi_pre[gpmi_pre==-9999.9] = np.nan
gpmi_pre = np.ma.masked_invalid(gpmi_pre)



print ('IMERG Rain Shape:',gpmi_pre.shape)
print ('IMERG Lat Shape:', gpmi_lat.shape)
print ('IMERG Lon Shape:', gpmi_lon.shape)
print ('IMERG Rain Range:',np.nanmin(gpmi_pre),' bis ',np.nanmax(gpmi_pre))

## Plot
## ----

gpmvmax = np.nanmax(gpmi_pre)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
pm2 = plt.pcolormesh(gpmi_lon, gpmi_lat, gpmi_pre,vmin=0,vmax=10)#,vmin=0,vmax=maxv)
cbar = plt.colorbar(pm2, shrink=0.75)
plot_borders(ax)
plt.xlabel("lon")
plt.ylabel("lat")
plt.xlim(-180,180)
plt.ylim(-90,90)
plt.grid(True)


plt.show()

