"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!

"""


"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

"""
import numpy as np
import matplotlib.pyplot as plt
import wradlib  
import wradlib as wrl
import pandas as pd
import datetime as dt
import pcc as pcc

from pcc import boxpol_pos
from pcc import plot_radar

bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']


############################################### Zeitstempel nach YYYYMMDDhhmmss
#Atime = pd.date_range('01/01/2014', periods=2, freq='5min')
#A = pd.Timestamp('20120501120500')


from pcc import zeitschleife as zt

zeit = zt(2017,1,5,2,25,0,
          2017,1,5,2,35,0)

print zeit

for ij in range(len(zeit)):
    print ij
    ZP = zeit[ij]
    print ZP
#ZP = '20150513180500'
    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]



    ################################################### Read RADOLAN GK Koordinaten

    iii = 0
    pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
            str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
            str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

    pfad_radolan = pfad[:-3]


    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)


    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

    radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
    radolan_zeit_sav = rwattrs['datetime'].strftime("%Y%m%d-%H%M%S")

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]

    Z = wradlib.trafo.idecibel(rwdata)

    # Marshall and Palmer 1948
    rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)

    ########################################################### Landgrenzenfunktion

    def plot_borders(ax):

        from osgeo import osr
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)
        india = osr.SpatialReference()
        # asia south albers equal area conic
        india.ImportFromEPSG(102028)

        proj_gk = osr.SpatialReference()
        proj_gk.ImportFromEPSG(31466)
        proj_ll = osr.SpatialReference()
        proj_ll.ImportFromEPSG(4326)
        gk3 = wradlib.georef.epsg_to_osr(31467)
        proj_stereo = wrl.georef.create_osr("dwd-radolan")
        proj_wgs = osr.SpatialReference()
        proj_wgs.ImportFromEPSG(4326)

        # country list
        countries = ['Germany']#,'France','Denmark', 'Netherlands', 'Poland']
        # open the input data source and get the layer
        filename = wradlib.util.get_wradlib_data_file('/automount/db01/python/data/NED/10m/cultural/10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp')
        dataset, inLayer = wradlib.io.open_shape(filename)
        # iterate over countries, filter accordingly, get coordinates and plot
        for item in countries:
            #print item
            # SQL-like selection syntax
            fattr = "(name='"+item+"')"
            inLayer.SetAttributeFilter(fattr)
            # get borders and names
            borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')

            for j in range(borders.shape[0]):
                bu = np.array(borders[j].shape)
                a = np.array(bu.shape)

                if a==1:
                    for i in range(0,borders[j].shape[0],1):
                        bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                        bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                        wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
                if a==2:    #richtig
                    bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                    bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                    wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
        ax.autoscale()


    dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/'
                                               'germany/vg250_0101.gk3.shape.ebenen'
                                               '/vg250_ebenen/vg250_l.shp')

    import matplotlib.cm as cm
    my_cmap = cm.get_cmap('jet',40)
    my_cmap.set_under('lightgrey')
    my_cmap.set_over('darkred')



    cmap2 = pcc.get_miub_cmap() #' bei Reflektivitat'
    ########################################################################## PLOT

    ff = 15
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(111, aspect='equal')
    plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    #plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
    cb = plt.colorbar(shrink=0.8)
    cb.set_label("Rainrate (mm/h)",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plot_borders(ax1)
    plt.title('RADOLAN Rainrate: \n'+ radolan_zeit + 'UTC',fontsize=ff)

    #plot_ocean(ax1)
    plt.xlabel("x [km] ",fontsize=ff)
    plt.ylabel("y [km]  ",fontsize=ff)
    #plt.xticks(fontsize=0)
    #plt.yticks(fontsize=0)
    plt.grid(color='r')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)
    plt.tight_layout()

    plot_radar(blon, blat, ax1, reproject=True)

    plt.savefig('/home/velibor/shkgpm/plot/radolan/rx_'+ radolan_zeit_sav+ '.png')
    plt.show()
    plt.close()
    #plt.show()
