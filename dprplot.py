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
ZP = '20141007023500'; gpm_time = '2014-10-07, 02:36 UTC'
#ZP = '20170113001000'; gpm_time = '2017-01-13, 00:12 UTC'
#ZP = '20140609'; gpm_time = '2014-06-09'
#ZP = '20140629'; gpm_time = '2014-06-29'
#ZP = '20140826'; gpm_time = '2014-08-26'
#ZP = '20160904'; gpm_time = '2016-09-04'
#ZP = '20140921'; gpm_time = '2014-09-21'
#ZP = '20141016'; gpm_time = '2014-10-16'
#ZP = '20150128'; gpm_time = '2015-01-28'
#ZP = '20150427'; gpm_time = '2015-04-27'
#ZP = '20160405'; gpm_time = '2016-04-05'
#ZP = '20160607'; gpm_time = '2016-06-07'




year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
ye = ZP[2:4]

## GPM DPR
## ----------------------------

pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]


gpmdprs = h5py.File(pfad_gprof_g, 'r')

# Ku Normalscan
ku_lat = np.array(gpmdprs['NS']['Latitude'])
ku_lon = np.array(gpmdprs['NS']['Longitude'])
ku_pp = np.array(gpmdprs['NS']['SLV']['precipRateNearSurface'])
ku_z = np.array(gpmdprs['NS']['SLV']['zFactorCorrectedNearSurface'])
ku_z[ku_z==-9999.9] = np.NaN
ku_pp[ku_pp==-9999.9] = np.NaN

# Ka Highintensive scan
ka_lat = np.array(gpmdprs['HS']['Latitude'])
ka_lon = np.array(gpmdprs['HS']['Longitude'])
ka_pp = np.array(gpmdprs['HS']['SLV']['precipRateNearSurface'])
ka_z = np.array(gpmdprs['HS']['SLV']['zFactorCorrectedNearSurface'])
ka_z[ka_z==-9999.9] = np.NaN
ka_pp[ka_pp==-9999.9] = np.NaN

# Ka+ku Matched scan
kaku_lat = np.array(gpmdprs['MS']['Latitude'])
kaku_lon = np.array(gpmdprs['MS']['Longitude'])
kaku_pp = np.array(gpmdprs['MS']['SLV']['precipRateNearSurface'])
kaku_z = np.array(gpmdprs['MS']['SLV']['zFactorCorrectedNearSurface'])
kaku_z[kaku_z==-9999.9] = np.NaN
kaku_pp[kaku_pp==-9999.9] = np.NaN



# CUT Ueber Deutschland
from pcc import cut_the_swath
#ku_lon, ku_lat, ku_pp = cut_the_swath(ku_lon,ku_lat,ku_pp)
#ka_lon, ka_lat, ka_pp = cut_the_swath(ka_lon,ka_lat,ka_pp)
#kaku_lon, kaku_lat, kaku_pp = cut_the_swath(kaku_lon,kaku_lat,kaku_pp)



proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

ku_x, ku_y = wradlib.georef.reproject(ku_lon, ku_lat, projection_target=proj_stereo , projection_source=proj_wgs)
ka_x, ka_y = wradlib.georef.reproject(ka_lon, ka_lat, projection_target=proj_stereo , projection_source=proj_wgs)
kaku_x, kaku_y = wradlib.georef.reproject(kaku_lon, kaku_lat, projection_target=proj_stereo , projection_source=proj_wgs)



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

ax1 = fig.add_subplot(2,3,1, aspect='equal')
pm1 = plt.pcolormesh(ku_x, ku_y,np.ma.masked_invalid(ku_pp),
                     cmap=my_cmap,
                     vmin=0.1,
                     vmax=10
                     )
plt.plot(ka_x[:,0],ka_y[:,0], color='red',lw=1)
plt.plot(ka_x[:,-1],ka_y[:,-1], color='red',lw=1)

plt.plot(kaku_x[:,0],kaku_y[:,0], color='red',lw=1)
plt.plot(kaku_x[:,-1],kaku_y[:,-1], color='red',lw=1)

cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ku: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax1)
plot_radar(boxlon, boxlat, ax1, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)

ax2 = fig.add_subplot(2,3,3, aspect='equal')
pm2 = plt.pcolormesh(ka_x, ka_y,np.ma.masked_invalid(ka_pp),
                     cmap=my_cmap,
                     vmin=0.1,
                     vmax=10
                     )
cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ka: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax2)
plot_radar(boxlon, boxlat, ax2, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)

ax3 = fig.add_subplot(2,3,2, aspect='equal')
pm3 = plt.pcolormesh(kaku_x, kaku_y,np.ma.masked_invalid(kaku_pp),
                     cmap=my_cmap,
                     vmin=0.01,
                     vmax=10,
                     )
cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ku + Ka: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax3)
plot_radar(boxlon, boxlat, ax3, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)


ax4 = fig.add_subplot(2,3,4, aspect='equal')
pm4 = plt.pcolormesh(ku_x, ku_y,np.ma.masked_invalid(ku_z),
                     cmap=my_cmap2,
                     vmin=0,
                     vmax=50

                     )
plt.plot(ka_x[:,0],ka_y[:,0], color='red',lw=1)
plt.plot(ka_x[:,-1],ka_y[:,-1], color='red',lw=1)

plt.plot(kaku_x[:,0],kaku_y[:,0], color='red',lw=1)
plt.plot(kaku_x[:,-1],kaku_y[:,-1], color='red',lw=1)

cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname2,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ku: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax4)
plot_radar(boxlon, boxlat, ax4, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)

ax5 = fig.add_subplot(2,3,6, aspect='equal')
pm5 = plt.pcolormesh(ka_x, ka_y,np.ma.masked_invalid(ka_z),
                     cmap=my_cmap2,
                     vmin=0,
                     vmax=50
                     )
cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname2,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ka: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax5)
plot_radar(boxlon, boxlat, ax5, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)

ax5 = fig.add_subplot(2,3,5, aspect='equal')
pm5 = plt.pcolormesh(kaku_x, kaku_y,np.ma.masked_invalid(kaku_z),
                     cmap=my_cmap2,
                     vmin=0,
                     vmax=50,
                     )
cb = plt.colorbar(shrink=0.5,extend='max')
cb.set_label(cbname2,fontsize=fft)
cb.ax.tick_params(labelsize=fft)
plt.title('GPM DPR Ku + Ka: \n'+ gpm_time ,fontsize=fft)
plot_borders(ax5)
plot_radar(boxlon, boxlat, ax5, reproject=True)
plt.grid(color='r')
plt.tight_layout()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
#plt.ylim(-4400,-4100)
#plt.xlim(-350,-80)

plt.show()



