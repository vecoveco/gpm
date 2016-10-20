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
import wradlib as wrl
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

pfad = ('/home/velibor/GPMDATA/arthurhou.pps.eosdis.nasa.gov/gpmdata/2014/06/09/imerg/3B*.HDF5')
pfad_imerg_raw = sorted(glob.glob(pfad))

for ii in range(0,len(pfad_imerg_raw),1):
    pfad_imerg = pfad_imerg_raw[ii]
    print (pfad_imerg[-39:-15])

###############################################
    year = pfad_imerg[-39:-35]
    m = pfad_imerg[-35:-33]
    d = pfad_imerg[-33:-31]
    ht = pfad_imerg[-29:-27]
    mt = pfad_imerg[-27:-25]
    st = pfad_imerg[-25:-23]

    ht_end = pfad_imerg[-21:-19]
    mt_end = pfad_imerg[-19:-17]
    st_end = pfad_imerg[-17:-15]

    print (year, m, d, ht, mt , st)
    ppi_datapath=('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",00.mvol")


    ppi=h5py.File(ppi_datapath,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

    ZH = data['SCAN0']['ZH']['data']
    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']
    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']
    Z = wradlib.trafo.idecibel(ZH)
    R = wradlib.zr.z2r(Z, a=200., b=1.6)
####################################################

    gpmi = h5py.File(pfad_imerg,'r')
    #gpmp = h5py.File(pfadg_imerg,'r')

    print ("GPM Data:",gpmi.keys())

    gpmi_lat = gpmi[u'Grid'][u'lat']
    gpmi_lon = gpmi[u'Grid'][u'lon']
    gpmi_pre = gpmi[u'Grid'][u'IRprecipitation']

    gpmi_pre = np.array(gpmi_pre).transpose()
    gpmi_pre[gpmi_pre==-9999.9] = np.nan
    gpmi_pre = np.ma.masked_invalid(gpmi_pre)



    ## Plot
    ## ----
    bonn_lat1 = 49.9400
    bonn_lat2 = 51.3500
    bonn_lon1 = 6.40000
    bonn_lon2 = 8.10000
    #lim1 [xmin, xmax, ymin,ymax]
    #limit= np.array([-140, -90, -10, 30])#pazifik
    #limit= np.array([65, 100, 5, 40]) #indien
    #limit= np.array([5, 9, 49, 53]) #bonn
    #limit= np.array([-140, -80, -5, 30])#pazifik
    #limit= np.array([2,18,44,58]) #deutschland
    #limit= np.array([-180,180,-90, 90]) #Welt
    limit= np.array([bonn_lon1,bonn_lon2,bonn_lat1,bonn_lat2])

    fig = plt.figure()
    plt.subplot(211)
    ax = fig.add_subplot(211, aspect='equal')
    pm2 = plt.pcolormesh(gpmi_lon, gpmi_lat, np.ma.masked_invalid(gpmi_pre),vmin=0,vmax=10)#,vmin=0,vmax=maxv)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("IMERG RainRate [mm/h]")
    plot_borders(ax)
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.xlim(limit[0], limit[1])
    plt.ylim(limit[2], limit[3])
    plt.title(str(pfad_imerg[-39:-15]))
    plt.grid(True)


    plt.subplot(212)  # ==== RainRate boxpol ==== #
    ax1, pm2 = wradlib.vis.plot_ppi(R,r,az,vmin=0,vmax=10)
    cbar = plt.colorbar(pm2, shrink=0.75)
    cbar.set_label("RainRate Boxpol [mm/h]")
    plt.xlim((-101000,101000))
    plt.ylim((-101000,101000))
    plt.xticks(())
    plt.yticks(())
    plt.xlabel("X Range ")
    plt.ylabel("Y Range ")
    plt.title(ppi_datapath[-28:-8])

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
    plt.savefig('/home/velibor/GPMPLOT/3B' + str(pfad_imerg[-39:-15]) + '.png')
    plt.close()
    print (pfad_imerg[-39:-15])



print ("FERTIG!")