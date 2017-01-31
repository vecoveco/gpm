"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
import wradlib as wrl
from osgeo import osr


#ZP = '20160805055000'; gpm_time = '2016-08-05 T: 054700 UTC'
#ZP = '20160607155500'; gpm_time = '2016-06-07 T: 155500 UTC'
#ZP = '20160405174500'; gpm_time = '2016-04-05 T: 174500 UTC'
#ZP = '20141007023500'; gpm_time = '2014-10-07, 02:36 UTC'
#ZP = '20170113001000'; gpm_time = '2017-01-13, 00:12 UTC'
#ZP = '20140609'; gpm_time = '2014-06-09'
#ZP = '20140629'; gpm_time = '2014-06-29'
#ZP = '20140826'; gpm_time = '2014-08-26'
#ZP = '20160904'; gpm_time = '2016-09-04'
ZP = '20140921'; gpm_time = '2014-09-21'
#ZP = '20141016'; gpm_time = '2014-10-16'
#ZP = '20150128'; gpm_time = '2015-01-28'
#ZP = '20150427'; gpm_time = '2015-04-27'
#ZP = '20160405'; gpm_time = '2016-04-05'
#ZP = '20160607'; gpm_time = '2016-06-07'




year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
ye = ZP[2:4]

## GPM DPR
## ----------------------------

pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]


gpmdprs = h5py.File(pfad_gprof_g, 'r')

# Ku Normalscan
gprof_lat = np.array(gpmdprs['S1']['Latitude'])
gprof_lon = np.array(gpmdprs['S1']['Longitude'])
gprof_pp = np.array(gpmdprs['S1']['surfacePrecipitation'])
gprof_pp[gprof_pp==-9999.9] = np.NaN




#Wippen
wippe_lat = np.array([50.7455,50.7294,50.7114,50.6886,50.7689,50.6911,50.6696,50.6630,
             50.6815,50.6570,50.7118,50.7369,50.7502,50.7517,50.7269,50.7238,
             50.6392,50.7049,50.7269,50.7092,50.7127,50.6773])
wippe_lon = np.array([7.07457,7.08032,7.12037,7.08439,7.06215,7.15079,7.18338,7.14733,
             7.13585,7.12403, 7.14917,7.12871,7.20207,7.16625,7.18750,
             7.14808,7.12054,7.17520,7.04964,7.07456,7.06453,7.06609])

# CUT Ueber Deutschland
from pcc import cut_the_swath
#ku_lon, ku_lat, ku_pp = cut_the_swath(ku_lon,ku_lat,ku_pp)
#ka_lon, ka_lat, ka_pp = cut_the_swath(ka_lon,ka_lat,ka_pp)
#kaku_lon, kaku_lat, kaku_pp = cut_the_swath(kaku_lon,kaku_lat,kaku_pp)



proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gprof_x, gprof_y = wradlib.georef.reproject(gprof_lon, gprof_lat, projection_target=proj_stereo , projection_source=proj_wgs)
w_x, w_y = wradlib.georef.reproject(wippe_lon, wippe_lat, projection_target=proj_stereo , projection_source=proj_wgs)




fig = plt.figure(figsize=(12,8))
fft=15
lll = np.arange(-30,60,1)

cbname = 'Rainrate (mm/h)'
cbname2 = 'Reflectivity (dBZ)'
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
boxlat, boxlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']

from pcc import plot_borders
from pcc import plot_radar
from pcc import get_my_cmap
from pcc import get_miub_cmap
my_cmap = get_my_cmap()
my_cmap2 = get_miub_cmap()





ax1 = fig.add_subplot(1,2,1, aspect='equal')
pm1 = plt.pcolormesh(gprof_x, gprof_y,np.ma.masked_invalid(gprof_pp),
                     cmap=my_cmap,
                     vmin=0.1,
                     vmax=10
                     )
plt.plot(gprof_x[:,0],gprof_y[:,0], color='red',lw=1)
plt.plot(gprof_x[:,-1],gprof_y[:,-1], color='red',lw=1)

plt.plot(gprof_x[:,0],gprof_y[:,0], color='red',lw=1)
plt.plot(gprof_x[:,-1],gprof_y[:,-1], color='red',lw=1)

cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ku: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax1)
plot_radar(boxlon, boxlat, ax1, reproject=True)
plt.scatter(w_x,w_y, color='red')
plt.scatter(gprof_x,gprof_y, color='black')
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax1 = fig.add_subplot(1,2,2, aspect='equal')
pm1 = plt.pcolormesh(gprof_lon, gprof_lat,np.ma.masked_invalid(gprof_pp),
                     cmap=my_cmap,
                     vmin=0.1,
                     vmax=10
                     )

plt.plot(gprof_lon[:,0],gprof_lat[:,0], color='red',lw=1)
plt.plot(gprof_lon[:,-1],gprof_lat[:,-1], color='red',lw=1)

plt.plot(gprof_lon[:,0],gprof_lat[:,0], color='red',lw=1)
plt.plot(gprof_lon[:,-1],gprof_lat[:,-1], color='red',lw=1)

cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.scatter(wippe_lon,wippe_lat, color='red')
plt.scatter(gprof_lon,gprof_lat, color='black')
#plt.xlim(7.0,7.3)
#plt.ylim(50.5,50.9)
plt.title('GPM DPR Ku: \n'+ gpm_time ,fontsize=fft)
#plot_borders(ax1)
#plot_radar(boxlon, boxlat, ax1, reproject=True)
plt.grid(color='r')
plt.tight_layout()

plt.show()





