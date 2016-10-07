"""Einlesen und darstellen von GPM IMERG Dateien"""

#wget -r --user=bregovic@gmx.de --password=bregovic@gmx.de ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/2015/01/01/imerg/

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

#FTP:::

pfad = ('/home/velibor/GPMDATA/arthurhou.pps.eosdis.nasa.gov/gpmdata/2015/01/01/imerg/3B*.HDF5')
pfad_imerg_raw = sorted(glob.glob(pfad))
print pfad_imerg_raw
#Todo: Schleife  Tage
for ii in range(0,len(pfad_imerg_raw),1):
    pfad_imerg = pfad_imerg_raw[ii]
    #pfadg = ('/user/velibor/SHKGPM/data/imerg/2A*.HDF5')
    #pfadg_imerg = sorted(glob.glob(pfadg))
    #pfadg_imerg = pfadg_imerg[0]

    gpmi = h5py.File(pfad_imerg,'r')
    #gpmp = h5py.File(pfadg_imerg,'r')

    print ("GPM Data:",gpmi.keys())

    gpmi_lat = gpmi[u'Grid'][u'lat']
    gpmi_lon = gpmi[u'Grid'][u'lon']
    #gpmi_pre = gpmi[u'Grid'][u'HQprecipitation']
    #gpmi_pre = gpmi[u'Grid'][u'precipitationUncal']
    #gpmi_pre = gpmi[u'Grid'][u'probabilityLiquidPrecipitation']
    gpmi_pre = gpmi[u'Grid'][u'IRprecipitation']


    #gpmp_lat = gpmp[u'S1'][u'Latitude']
    #gpmp_lon = gpmp[u'S1'][u'Longitude']
    ##gpmi_pre = gpmi[u'Grid'][u'HQprecipitation']
    ##gpmi_pre = gpmi[u'Grid'][u'precipitationUncal']
    ##gpmi_pre = gpmi[u'Grid'][u'probabilityLiquidPrecipitation']
    #gpmp_pre = gpmp[u'S1'][u'surfacePrecipitation']

    #gpmp_pre = np.array(gpmp_pre)#.transpose()
    #gpmp_pre[gpmp_pre==-9999.9] = np.nan
    #gpmp_pre = np.ma.masked_invalid(gpmp_pre)

    gpmi_pre = np.array(gpmi_pre).transpose()
    gpmi_pre[gpmi_pre==-9999.9] = np.nan
    gpmi_pre = np.ma.masked_invalid(gpmi_pre)



    ## Plot
    ## ----

    #lim1 [xmin, xmax, ymin,ymax]
    #limit= np.array([-140, -90, -10, 30])#pazifik
    #limit= np.array([65, 100, 5, 40]) #indien
    #limit= np.array([5, 9, 49, 53]) #bonn
    #limit= np.array([-140, -80, -5, 30])#pazifik
    limit= np.array([2,18,44,58]) #deutschland
    #limit= np.array([-180,180,-90, 90]) #Welt

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    pm2 = plt.pcolormesh(gpmi_lon, gpmi_lat, np.ma.masked_invalid(gpmi_pre),vmin=0,vmax=10)#,vmin=0,vmax=maxv)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("IMERG RainRate [mm/h]")
    plot_borders(ax)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.xlim(limit[0], limit[1])
    plt.ylim(limit[2], limit[3])
    plt.title(str(pfad_imerg[-39:-14]))
    plt.grid(True)

    #ax = fig.add_subplot(212, aspect='equal')
    #pm2 = plt.pcolormesh(gpmp_lon, gpmp_lat, np.ma.masked_invalid(gpmp_pre),vmin=0,vmax=10)#,vmin=0,vmax=maxv)
    #cbar = plt.colorbar(pm2, shrink=0.75)
    #plot_borders(ax)
    #plt.xlabel("lon")
    #plt.ylabel("lat")
    #plt.xlim(limit[0], limit[1])
    #plt.ylim(limit[2], limit[3])
    #plt.title(str(pfadg_imerg[-60:-1]))
    #plt.grid(True)
    #plt.tight_layout()
    #plt.show()

    plt.tight_layout()
    plt.savefig('/home/velibor/GPMPLOT/3B' + str(pfad_imerg[-39:-14]) + '.png')
    plt.close()
    print (pfad_imerg[-39:-14])


print ("FERTIG!")