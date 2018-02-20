"""


Plott DPR VCUT etc...


"""

"""
Analyse of GPM DPR BB

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import wradlib
import wradlib as wrl
from osgeo import osr
import h5py
import glob
from pcc import get_my_cmap
from pcc import get_miub_cmap
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']


from pcc import plot_borders
from pcc import plot_radar


#dates ='20140729'; jj = 1; cut1, cut2 = 24, 3254
#dates ='20151216'; jj = 0; cut1, cut2 = 41, 4660
#dates ='20140806'; jj = 1; cut1, cut2 = 41, 4660
#dates ='20170211'; jj = 0; cut1, cut2 = 17, 4660
#dates ='20161109'; jj = 1; cut1, cut2 = 25, 4660
dates ='20160209'; jj = 1; cut1, cut2 = 25, 4650

#def gpm_bb(dates, pn=0):
zt=dates

pfad = ('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.'+zt+'*.HDF5')
dpr_pfad = sorted(glob.glob(pfad))[jj]

print dpr_pfad

scan  = 'NS'#NS' #or MS

# Einlesen
dpr = h5py.File(dpr_pfad, 'r')
dpr_lat=np.array(dpr[scan]['Latitude'])
dpr_lon=np.array(dpr[scan]['Longitude'])
dpr_pp=np.array(dpr[scan]['SLV']['zFactorCorrected'])
dpr_pp[dpr_pp<0]= np.nan

dpr_pp_surf=np.array(dpr[scan]['SLV']['zFactorCorrectedNearSurface'])
dpr_pp_surf[dpr_pp_surf<0]= np.nan


dpr_bbh=np.array(dpr[scan]['CSF']['heightBB'], dtype=float)
dpr_bbh[dpr_bbh<0]= np.nan

dpr_bbh_flag = np.array(dpr[scan]['CSF']['flagBB'], dtype=float)
dpr_bbh_flag[dpr_bbh_flag<0]= np.nan

dpr_bbh_qual = np.array(dpr[scan]['CSF']['qualityBB'], dtype=float)
dpr_bbh_qual[dpr_bbh_qual<0]= np.nan


dpr_bbw=np.array(dpr[scan]['CSF']['widthBB'], dtype=float)
dpr_bbw[dpr_bbw<0]= np.nan

dpr_time = dpr['NS']['ScanTime']

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
blat, blon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']

dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat, projection_target=proj_stereo , projection_source=proj_wgs)
bonnlon, bonnlat = wradlib.georef.reproject(bonnlon, bonnlat, projection_target=proj_stereo , projection_source=proj_wgs)


lon0, lat0, radius = bonnlon, bonnlat, 100
r = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = r < radius

lat = dpr_lat[position]
lon = dpr_lon[position]

dpr_pp[np.where(r > radius)]=np.nan
pp=dpr_pp

dpr_pp_surf[np.where(r > radius)]=np.nan

dpr_bbw[np.where(r > radius)]=np.nan
dpr_bbh[np.where(r > radius)]=np.nan



h = np.arange(150,4800,150)
if scan=='HS':
    hdpr = 1000 * (np.arange(88,0,-1)*0.250)

else:
    hdpr = 1000 * (np.arange(176,0,-1)*0.125)

hhh = np.array(pp.shape[0]*pp.shape[1]*list(hdpr))
ppp = pp.reshape(pp.shape[0]*pp.shape[1]*pp.shape[2])

maske = ~np.isnan(hhh) & ~np.isnan(ppp)


#cut1, cut2 = 24, 3254

fig = plt.figure(figsize=(14,12))
#zzz = str(jahr)+'-'+str(monat)+'-'+str(tag)+'--'+str(stunde)+':'+str(minute)+' UTC'
#fig.suptitle(zzz + ' UTC')

###################
ax1 = fig.add_subplot(221, aspect='auto')
#plt.subplot(2,2,1)
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp_surf), vmin=0, vmax=50, cmap=get_miub_cmap())
cbar = plt.colorbar()
cbar.set_label('Ref. in dbz')
plot_borders(ax1)
plot_radar(blon, blat, ax1, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')
plt.plot(dpr_lon[:,cut1],dpr_lat[:,cut1], color='red',lw=2,ls='--')
plt.plot(dpr_lon[cut2,:],dpr_lat[cut2,:], color='green',lw=2,ls='--')


ax1 = plt.scatter(bonnlon, bonnlat, c=50 ,s=50, color='red')


plt.grid()
plt.xlim(-350,-100)
plt.ylim( -4350,-4100)

##################exit()
hhh = hhh/1000.
ax2 = fig.add_subplot(222, aspect='auto')
ax2.hist2d(ppp[maske],hhh[maske], bins=30, cmap=get_my_cmap(), vmin=0.1)

#plt.plot(np.nanmax(pp[:,:],axis=0),hdpr, color='red', lw=2)
ax2.plot(np.nanmean(pp[:,:,:],axis=(0,1)),hdpr/1000., color='red', lw=2)
plt.plot(np.nanmedian(pp[:,:,:],axis=(0,1)),hdpr/1000., color='green', lw=2)

cbar = plt.colorbar()
cbar.set_label('number of samples')

#plt.title('DPR Ref. in Box')
plt.xlabel('Reflectivity  (dBZ)')
plt.ylabel('Height (km)')
plt.grid()
plt.xticks()
plt.yticks()

#plt.ylim(0,6000)
#plt.xlim(0,50)
##################
#print np.uniforn(bbh)
#mini = np.nanmin(bbh[bbh>0])
ax3 = fig.add_subplot(223, aspect='auto')
#dpr_lon[:,cut], hdpr,
plt.pcolormesh(dpr_lon[:,cut1],hdpr,np.ma.masked_invalid(dpr_pp[:,cut1,:]).T , vmin=0, vmax=50, cmap=get_miub_cmap())#np.nanmin(dpr_bbh[dpr_bbh>0])
cbar = plt.colorbar()
cbar.set_label('-')
plt.title('vcut', color='red')


plt.grid()
#plt.title('BB Hight')
plt.xlim(-350,-100)
plt.ylim(0, 8000)

##################
ax4 = fig.add_subplot(224, aspect='auto')
#dpr_lon[:,cut], hdpr,
plt.pcolormesh(dpr_lon[cut2,:],hdpr,np.ma.masked_invalid(dpr_pp[cut2,:,:]).T, vmin=0, vmax=50, cmap=get_miub_cmap())#np.nanmin(dpr_bbh[dpr_bbh>0])
cbar = plt.colorbar()
cbar.set_label('-')
plt.title('vcut', color='green')
plt.xlim(-350,-100)
plt.ylim(0, 8000)

plt.grid()

plt.tight_layout()
plt.show()
#plt.savefig('/automount/ags/velibor/plot/BB/_'+scan+'dprbb_'+str(dates)+'.png' )
#plt.close()

