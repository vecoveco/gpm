"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats
import wradlib as wrl
from osgeo import osr



ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
TH_rain= 0.2

# Zeitstempel nach YYYYMMDDhhmmss
#ZP = '20141007023500'#'20141007023500'#'20161024232500'#'20150427223500' #'20141007023500'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'
#ZP = '20170203005500'
ZP = '20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------

pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]
rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5






pfadz = ('/automount/radar/dwd/rz/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rz_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolanz = pfadz
print pfad_radolanz


rw_filenamez = wradlib.util.get_wradlib_data_file(pfad_radolanz)
rzdata, rzattrs = wradlib.io.read_RADOLAN_composite(rw_filenamez)




radolan_grid_xyz = wradlib.georef.get_radolan_grid(900,900)
xz = radolan_grid_xyz[:,:,0]
yz = radolan_grid_xyz[:,:,1]

rzdata = np.ma.masked_equal(rzdata, -9999)

Zr = wradlib.trafo.idecibel(rwdata)
rrr = wradlib.zr.z2r(Zr, a=200., b=1.6)

ff = 15

plt.subplot(1,2,1)
plt.pcolormesh(x,y,rrr, vmin=0, vmax=10)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
plt.subplot(1,2,2)
plt.pcolormesh(xz,yz,rzdata, vmin=0, vmax=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)

plt.show()






