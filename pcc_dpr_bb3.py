"""
Analyse of GPM DPR BB

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import wradlib
import satlib as sl
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


dates ='20141007'

zt=dates

pfad = ('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.'+zt+'*.HDF5')
dpr_pfad = sorted(glob.glob(pfad))[0]

print dpr_pfad

scan  = 'MS' #or MS

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
dpr_bbw=np.array(dpr[scan]['CSF']['widthBB'], dtype=float)
dpr_bbw[dpr_bbw<0]= np.nan



proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat, projection_target=proj_stereo , projection_source=proj_wgs)
bonnlon, bonnlat = wradlib.georef.reproject(bonnlon, bonnlat, projection_target=proj_stereo , projection_source=proj_wgs)

print '-------->', bonnlon, bonnlat


lon0, lat0, radius = bonnlon, bonnlat, 100
r = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
position = r < radius
#plt.scatter(dpr_lon[position], dpr_lat[position])
#plt.scatter(lon0,lat0,c=40, s=40, color='red')
#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
#plt.grid()
#plt.show()


#position = np.where((dpr_lat<l1) & (dpr_lat>l2) & (dpr_lon < k1) & (dpr_lon > k2))


lat = dpr_lat[position]
lon = dpr_lon[position]
#pp = dpr_pp[position]

dpr_pp[np.where(r > radius)]=np.nan
pp=dpr_pp
print pp.shape

dpr_pp_surf[np.where(r > radius)]=np.nan

dpr_bbw[np.where(r > radius)]=np.nan
dpr_bbh[np.where(r > radius)]=np.nan

#bbw = dpr_bbw[position]
#bbh = dpr_bbh[position]

#l2, l1 = -190, -250
#k2, k1 = -4210, -4270
#pos = np.where((dpr_lat<l1) & (dpr_lat>l2) & (dpr_lon < k1) & (dpr_lon > k2))
#dpr_time = dpr['NS']['ScanTime']

#stunde = np.array(dpr_time['Hour'])[pos[0]][0]
#minute = np.array(dpr_time['Minute'])[pos[0]][0]
#sekunde = np.array(dpr_time['Second'])[pos[0]][0]

#jahr = np.array(dpr_time['Year'])[pos[0]][0]
#monat = np.array(dpr_time['Month'])[pos[0]][0]
#tag = np.array(dpr_time['DayOfMonth'])[pos[0]][0]

#zeit =  (str(jahr)+'.'+str(monat)+'.'+str(tag) + ' -- ' + str(stunde)+':'+str(minute)+':'+str(sekunde))
#print zeit

h = np.arange(150,4800,150)
if scan=='HS':
    hdpr = 1000 * (np.arange(88,0,-1)*0.250)

else:
    hdpr = 1000 * (np.arange(176,0,-1)*0.125)

pp.shape
hhh = np.array(pp.shape[0]*pp.shape[1]*list(hdpr))
ppp = pp.reshape(pp.shape[0]*pp.shape[1]*pp.shape[2])

maske = ~np.isnan(hhh) & ~np.isnan(ppp)




fig = plt.figure(figsize=(12,10))
#zzz = str(jahr)+'-'+str(monat)+'-'+str(tag)+'--'+str(stunde)+':'+str(minute)+' UTC'
#fig.suptitle(zzz + ' UTC')

###################
ax1 = fig.add_subplot(221, aspect='auto')
#plt.subplot(2,2,1)
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp_surf), vmin=np.nanmin(dpr_pp_surf), vmax=np.nanmax(dpr_pp_surf), cmap=get_miub_cmap())
cbar = plt.colorbar()
cbar.set_label('Ref. in dbz')
plot_borders(ax1)
plot_radar(blon, blat, ax1, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')


ax1 = plt.scatter(bonnlon, bonnlat, c=50 ,s=50, color='red')


plt.grid()
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

##################
ax2 = fig.add_subplot(222, aspect='auto')
plt.hist2d(ppp[maske],hhh[maske], bins=30, cmap=get_my_cmap(), vmin=0.1)
print pp.shape

#plt.plot(np.nanmax(pp[:,:],axis=0),hdpr, color='red', lw=2)
plt.plot(np.nanmean(pp[:,:,:],axis=(0,1)),hdpr, color='red', lw=2)
plt.plot(np.nanmedian(pp[:,:,:],axis=(0,1)),hdpr, color='green', lw=2)
cbar = plt.colorbar()
cbar.set_label('#')


plt.title('DPR Ref. in Box')
plt.xlabel('Reflectivity in dBZ')
plt.grid()
plt.xticks()
plt.yticks()

#plt.ylim(0,6000)
#plt.xlim(0,50)
##################
#print np.uniforn(bbh)
#mini = np.nanmin(bbh[bbh>0])

ax3 = fig.add_subplot(223, aspect='auto')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbh), vmin=np.nanmin(dpr_bbh[dpr_bbh>0]), vmax=np.nanmax(dpr_bbh), cmap='jet')
cbar = plt.colorbar()
cbar.set_label('Hight in m')

plot_borders(ax3)
plot_radar(blon, blat, ax3, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

ax1 = plt.scatter(bonnlon, bonnlat, c=50 ,s=50, color='red')
plt.grid()
plt.title('BB Hight')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

##################
ax4 = fig.add_subplot(224, aspect='auto')
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbw), vmin=np.nanmin(dpr_bbw[dpr_bbh>0]), vmax=np.nanmax(dpr_bbw), cmap='jet')
cbar = plt.colorbar()
cbar.set_label('Width in m')

plot_borders(ax4)
plot_radar(blon, blat, ax4, reproject=True, cband=False,col='black')
plt.plot(dpr_lon[:,0],dpr_lat[:,0], color='black',lw=1)
plt.plot(dpr_lon[:,-1],dpr_lat[:,-1], color='black',lw=1)
plt.plot(dpr_lon[:,dpr_lon.shape[1]/2],dpr_lat[:,dpr_lon.shape[1]/2], color='black',lw=1, ls='--')

ax1 = plt.scatter(bonnlon, bonnlat, c=50 ,s=50, color='red')
plt.grid()
plt.title('BB Width')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

plt.tight_layout()
plt.show()
#plt.savefig('/automount/ags/velibor/plot/BB/'+scan+'dprbb_'+str(zzz)+'.png' )
#plt.close()


#dates = ['20150128','20150330', '20150404','20151208','20151216','20160107','20160612','20161019','20161222',
#         '20141007','20140708', '20151015','20160209', '20160915','20161121', '20141008','20160601','20161024','20161109', '20140729']



    #pfad = ('/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.*.'+zt+'*.HDF5')

