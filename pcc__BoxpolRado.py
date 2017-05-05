"""

Vergleich Radolan und BOxPOl


"""

from satlib import read_rado
import matplotlib.pyplot as plt
import wradlib
import h5py
from osgeo import osr

zz = '20141007023000'

x, y, rwdata, rn = read_rado(zz, r_pro='rx')



def read_boxpol(zeitstempel):
    ZP = zeitstempel


    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]

    if year < '2015':
        print 'archive'
        sc = 'radar-archiv'
    else:
        sc = 'radar'

    pfad = ('/automount/'+sc+'/scans/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
            str(year)+'-'+str(m)+'-'+str(d)+'/ppi_1p5deg/'+str(year)+'-'+str(m)+'-'+
            str(d)+'--'+str(ht)+':'+str(mt)+':00,00.mvol')

    boxpol_filename = wradlib.util.get_wradlib_data_file(pfad)

    ppi = h5py.File(boxpol_filename,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(boxpol_filename)


    zh = data['SCAN0'][u'ZH']['data']
    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']

    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']

    return zh, r, az, lon_ppi, lat_ppi, alt_ppi

zh, r, az, lon_ppi, lat_ppi, alt_ppi = read_boxpol(zz)

"""proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
"""

import numpy as np


gky = -4235.233235191105,
gkx = -216.64772430049572







###############################################################################
radar_location = (lon_ppi, lat_ppi, alt_ppi) # (lon, lat, alt) in decimal degree and meters
elevation = 1.5 # in degree
azimuths = az # in degrees
ranges = r # in meters
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)


#ae = wradlib.georef.create_osr("aeqd", lon_0=radar_location[0], lat_0=radar_location[1])
#x, y = wradlib.georef.reproject(lon, lat, projection_target=ae)
proj_stereo = wradlib.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

bx, by = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((bx.ravel(), by.ravel())).transpose()


###############################################################################
grid_xy = np.vstack((x.ravel(), y.ravel())).transpose() # Radolangitter
xy=np.concatenate([bx.ravel()[:,None],by.ravel()[:,None]], axis=1)

gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([bx.mean(), by.mean()]), zh.ravel(), wradlib.ipol.Idw)#Linear, Idw, Nearest
gridded = np.ma.masked_invalid(gridded).reshape(x.shape)


from pcc import get_miub_cmap
from pcc import plot_radar


plt.subplot(2,2,3)
#wradlib.vis.plot_ppi(zh,r,az, vmin=0, vmax=50)
plt.pcolormesh(bx, by,zh, cmap=get_miub_cmap(), vmin=0, vmax=50)
plt.scatter(gkx, gky, s=25, color='black')

plt.subplot(2,2,2)
plt.pcolormesh(x,y,rwdata, vmin=0, vmax=50, cmap=get_miub_cmap())
plt.scatter(gkx, gky, s=25, color='black')

plt.xlim(-350, -100)
plt.ylim(-4400, -4100)

plt.subplot(2,2,1)
plt.pcolormesh(x,y,gridded, vmin=0, vmax=50, cmap=get_miub_cmap())
plt.scatter(gkx, gky, s=25, color='black')

plt.xlim(-350, -100)
plt.ylim(-4400, -4100)


plt.subplot(2,2,4)
a = gridded.copy()
b = rwdata.copy()
a[a<0]=np.nan
b[b<0]=np.nan
plt.scatter(a,b)
plt.show()