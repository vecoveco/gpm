"""Einlesen und darstellen von GPM IMERG Dateien"""

#wget -r --user=bregovic@gmx.de --password=bregovic@gmx.de ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/2015/01/01/imerg/

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import wradlib as wrl
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

## Einlesen von Parametern
## -----------------------------

bonn_lat1 = 49.9400
bonn_lat2 = 51.3500
bonn_lon1 = 6.40000
bonn_lon2 = 8.10000

limit= np.array([bonn_lon1,bonn_lon2,bonn_lat1,bonn_lat2])

my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')

## Einlesen von Radardaten
## -----------------------------
deg = ['ppi_1p5','ppi_3p4']

pfad = ('/automount/radar-archiv/scans/2014/2014-06/2014-06-09/' + str(deg[0]) +'deg/*')
pfadup = ('/automount/radar-archiv/scans/2014/2014-06/2014-06-09/' + str(deg[1]) +'deg/*')

ppi_datapath = sorted(glob.glob(pfad))
ppi_datapathup = sorted(glob.glob(pfadup))

for ii in range(0,len(ppi_datapath),1):
    pfad_ppi = ppi_datapath[ii]
    pfad_ppiup = ppi_datapathup[ii]
    print (pfad_ppi)

    rad_deg = [pfad_ppi, pfad_ppiup]
    subsub = [211,222]

    for jj in range(0,len(rad_deg),1):
        data, attrs = wradlib.io.read_GAMIC_hdf5(rad_deg[jj])

        ZH = data['SCAN0']['ZH']['data']
        r = attrs['SCAN0']['r']
        az = attrs['SCAN0']['az']
        lon_ppi = attrs['VOL']['Longitude']
        lat_ppi = attrs['VOL']['Latitude']
        alt_ppi = attrs['VOL']['Height']
        Z = wradlib.trafo.idecibel(ZH)
        R = wradlib.zr.z2r(Z, a=200., b=1.6)

        ## Plot

        fig = plt.figure()
        print jj
        plt.subplot(subsub[jj])
        ax1, pm2 = wradlib.vis.plot_ppi(R,r/1000,az,vmin=0.1,vmax=20, cmap=my_cmap)
        cbar = plt.colorbar(pm2, shrink=0.80)
        cbar.set_label("RainRate Boxpol [mm/h]")
        #plt.xticks(())
        #plt.yticks(())
        plt.xlim(-100,100)
        plt.ylim(-100,100)
        plt.xlabel("X-Range (km) ")
        plt.ylabel("Y-Range (km) ")
        plt.title('Rainrate Bonn: \n Time: ' + str(pfad_ppi[65:85])+' UTC')
        #fig.patch.set_visible(False)
        #ax1.axis('off')
        plt.hold(True)
    #plt.tight_layout()
    #plt.savefig('/home/velibor/GPMPLOT/Radar' + str(pfad_ppi[65:85])  +'.png')
    plt.close()
    print (pfad_ppi[65:85])

print ("FERTIG!")
