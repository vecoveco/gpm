
'''Dieses Program soll dazu dienen die
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Reflektivitat zu validieren.
Hier werden mehrere Ueberflug analysiert'''

#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import wradlib as wrl
import wradlib
import glob
import h5py


#data,attributes=wrl.io.read_GAMIC_hdf5('/automount/radar-archiv/scans/2014/2014-10/2014-10-07/ppi_1p5deg/2014-10-07--02:37:44,00.mvol')
ZP = '20141007023744'


year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]


try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]
    #ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/n_ppi_280deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

print ppi_datapath


# PPI BoxPol Daten einlesen
#---------------------------

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH = data['SCAN0']['ZH']['data']
PHIDP = data['SCAN0']['PHIDP']['data']
RHOHV = data['SCAN0']['RHOHV']['data']

r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']
th = attrs['SCAN0']['el']

A = ZH




r = np.arange(0, A.shape[1])
az = np.arange(0, A.shape[0])
# mask data array for better presentation
mask_ind = np.where(A <= np.nanmin(A))
#data[mask_ind] = np.nan
ma = np.ma.array(A, mask=np.isnan(A))

from pcc import get_miub_cmap as my_cmap

cgax, caax, paax, pm = wradlib.vis.plot_cg_ppi(ma, A, az, autoext=True,

                                           refrac=True, cmap=my_cmap(),vmin=0,vmax=50)

t = plt.title('BoXPol PPI')
t.set_y(1.05)
cbar = plt.gcf().colorbar(pm, pad=0.075, orientation='horizontal',shrink=0.7, extend='both')
#plt.text(1.0, 1.05, 'azimuth', transform=caax.transAxes, va='bottom',
#    ha='right')
cbar.set_label('Reflectivity (dBZ)')
plt.tick_params(
    axis='off',
    which='off',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off',fontsize=0)
plt.tight_layout()

plt.show()
'''
from pcc import get_miub_cmap as my_cmap
from pcc import get_my_cmap as my_cmap


import matplotlib.pyplot as plt
import numpy as np
# well, it's a wradlib example
import wradlib
from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
# reading in data, range and theta arrays from special rhi hdf5 file
file = '/automount/radar-archiv/scans/2014/2014-10/2014-10-07/n_rhi_lacros/2014-10-07--02:37:44,00.mvol'
#data, meta = wradlib.io.from_hdf5(file, dataset='data')
#r, meta = wradlib.io.from_hdf5(file, dataset='range')
#th, meta = wradlib.io.from_hdf5(file, dataset='theta')


boxpol_filename = wradlib.util.get_wradlib_data_file(file)
ppi=h5py.File(boxpol_filename,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(boxpol_filename)

print data[u'SCAN0'].keys()
print attrs['VOL'].keys()
print attrs['SCAN0'].keys()


zh = data['SCAN0'][u'ZH']['data']
phidp = data['SCAN0'][u'PHIDP']['data']
rhohv = data['SCAN0'][u'RHOHV']['data']
zv = data['SCAN0'][u'ZV']['data']
zdr = data['SCAN0'][u'ZDR']['data']
kdp = data['SCAN0'][u'KDP']['data']

r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
th = attrs['SCAN0']['el']

lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']


# mask data array for better presentation
#mask_ind = np.where(data <= np.nanmin(data))
#data[mask_ind] = np.nan
#ma = np.ma.array(data, mask=np.isnan(data))

# the simplest call, plot cg rhi in new window

cgax, caax, paax, pm = wradlib.vis.plot_cg_rhi(zdr, r=r, th=th, rf=1e3, refrac=False,
                                       subplot=111, cmap=my_cmap(), vmin=-0.1, vmax=2)
t = plt.title('BoXPol RHI')
t.set_y(1.05)
cgax.set_ylim(0,12)
cbar = plt.gcf().colorbar(pm, pad=0.05)
cbar.set_label('reflectivity [dBZ]')
#caax.set_xlabel('x_range [km]')
#caax.set_ylabel('y_range [km]')
#plt.text(1.0, 1.05, 'azimuth', transform=caax.transAxes, va='bottom',
#    ha='right')
gh = cgax.get_grid_helper()
# set theta to some nice values
#gh.grid_finder.grid_locator1 = FixedLocator([i for i in np.arange(0, 359, 5)])
locs = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.,
                15., 16., 17., 18., 20., 22., 25., 30., 35.,  40., 50., 60., 70., 80., 90.]
gh.grid_finder.grid_locator1 = FixedLocator(locs)
gh.grid_finder.tick_formatter1 = DictFormatter(dict([(i, r"${0:.0f}^\circ$".format(i)) for i in locs]))
plt.tight_layout()
plt.plot()
plt.show()
'''
'''
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
# ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')
from pcc import get_miub_cmap as my_cmap
from pcc import plot_radar
from pcc import boxpol_pos
import wradlib as wrl
from osgeo import osr

Pos = boxpol_pos()
bblon, bblat = Pos['lon_ppi'], Pos['lat_ppi']
gkx, gky = Pos['gkx_ppi'], Pos['gky_ppi']

# Pfad mit String
# ---------------

# Hohe von DPR
TH = 12 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


ZP = '20141007023744'
#ZP = '20140704134500'
#ZP = '20150225163500'
#ZP = '20150816070500'
#good#ZP = '20151208213500'
#good#ZP = '20151216024501'
#ZP = '20160503024500'

year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]



pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
print pfad_radar
#pfad_radar = sorted(glob.glob(pfad))
#print pfad_radar
pfad_radar_Ku = pfad_radar[0]

try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]
except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

print ppi_datapath


# PPI BoxPol Daten einlesen
#---------------------------

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH = data['SCAN0']['ZH']['data']
PHIDP = data['SCAN0']['PHIDP']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']

R = ZH
from wradlib.trafo import idecibel, decibel
R = idecibel(R)

#R = R + pia
#R[151:165]=np.nan



# DPR Einlesen
# ------------

gpmku = h5py.File(pfad_radar_Ku, 'r')
gpmku_HS=gpmku['NS']['SLV']
ku_lat=np.array(gpmku['NS']['Latitude'])			#(7934, 24)
ku_lon=np.array(gpmku['NS']['Longitude'])			#(7934, 24)
ku_pp=np.array(gpmku_HS['zFactorCorrectedNearSurface'])


wippe_lat = np.array([50.7455,50.7294,50.7114,50.6886,50.7689,50.6911,50.6696,50.6630,
             50.6815,50.6570,50.7118,50.7369,50.7502,50.7517,50.7269,50.7238,
             50.6392,50.7049,50.7269,50.7092,50.7127,50.6773])
wippe_lon = np.array([7.07457,7.08032,7.12037,7.08439,7.06215,7.15079,7.18338,7.14733,
             7.13585,7.12403, 7.14917,7.12871,7.20207,7.16625,7.18750,
             7.14808,7.12054,7.17520,7.04964,7.07456,7.06453,7.06609])

# Lon Lat Bestimmung
# ------------------
radars = [ku_pp]
rad_lat = [ku_lat]
rad_lon = [ku_lon]
rad_name = ['DPR']

ii = 0

dpr_pp = radars[ii]
dpr_lat = rad_lat[ii]
dpr_lon = rad_lon[ii]
radarname = rad_name[ii]

bonn_lat1 = 49.9400
bonn_lat2 = 51.3500
bonn_lon1 = 6.40000
bonn_lon2 = 8.10000

ilat= np.where((dpr_lat>49.9400) & (dpr_lat<51.3500))
ilon= np.where((dpr_lon>6.40000) & (dpr_lon<8.10000))
lonstart = ilon[0][0]
lonend = ilon[0][-1]
latstart = ilat[0][0]
latend = ilat[0][-1]
dpr_pp[dpr_pp==-9999] = np.nan
#dpr_pp = dpr_pp[:,:,hi]  # Untersete Schicht

radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

gk3 = wradlib.georef.epsg_to_osr(31467)
x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
xgrid, ygrid = wradlib.georef.reproject(dpr_lon[latstart:latend], dpr_lat[latstart:latend], projection_target=gk3)

grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=40,p=2)
gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)


gridded = decibel(gridded)

R = decibel(R)

# ON RADOLAN GRID
proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

box_x, box_y = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)
gpm_x, gpm_y = wradlib.georef.reproject(dpr_lon[latstart:latend], dpr_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)
wippe_x, wippe_y = wradlib.georef.reproject(wippe_lon, wippe_lat, projection_target=proj_stereo , projection_source=proj_wgs)

#Todo IM PLOT VERWENDEN !!!!!!!

# Plot
# ----
print '---------------------------------------------'
print wippe_x, wippe_y



################################################################Swap!
#rrr, ggg = ggg, rrr

ff = 20
cc = 0.5
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111, aspect='equal')#------------------------------------
#ax1, pm1 = wradlib.vis.plot_ppi(R,r,az,vmin=0.01,vmax=50, cmap=my_cmap())
pm1 = plt.pcolormesh(box_x, box_y, R, vmin=0, vmax=50, cmap=my_cmap())

plt.scatter(wippe_x, wippe_y, label='BoXPol Gauge Rain Tipper', s=200, marker='v')

#cb = plt.colorbar(pm1,shrink=cc)
#cb.set_label("Reflectivity (dBZ)",fontsize=ff)
#cb.ax.tick_params(labelsize=ff)

plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')

from pcc import plot_borders
plot_borders(ax1)
plot_radar(bblon, bblat, ax1, reproject=True, cband=False,col='black')

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')

plt.title('BoXPol Gauge Rain Tipper',fontsize=ff)

plt.grid(color='r')
plt.xlim(-321, -110)
plt.ylim(-4343,-4130)

plt.legend(loc='lower right', numpoints=1)

plt.tight_layout()
#plt.savefig('/home/velibor/shkgpm/plot/gpm_dpr_radolan_v2_'+ZP + '.png' )
#plt.close()
plt.show()
'''




