"""

Reading MRR Data

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import satlib as sl
import h5py



# DPR
dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20141007-S015721-E032951.003445.V04A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20160107-S120629-E133900.010562.V04A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170306-S151005-E164237.017160.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170223-S182621-E195854.016991.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170222-S113453-E130727.016971.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170401-S003803-E021037.017555.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170603-S235100-E012332.018550.V05A.HDF5'

h = np.arange(150,4800,150)

scan  = 'NS' #or MS

dpr = h5py.File(dpr_pfad, 'r')
dpr_lat=np.array(dpr[scan]['Latitude'])
dpr_lon=np.array(dpr[scan]['Longitude'])
dpr_pp=np.array(dpr[scan]['SLV']['zFactorCorrected'])
dpr_pp[dpr_pp<0]= np.nan

dpr_pp_surf=np.array(dpr[scan]['SLV']['zFactorCorrectedNearSurface'])
dpr_pp_surf[dpr_pp_surf<0]= np.nan

dpr_bbb=np.array(dpr[scan]['CSF']['binBBBottom'], dtype=float)
dpr_bbb[dpr_bbb<0]= np.nan
dpr_bbp=np.array(dpr[scan]['CSF']['binBBPeak'], dtype=float)
dpr_bbp[dpr_bbp<0]= np.nan
dpr_bbt=np.array(dpr[scan]['CSF']['binBBTop'], dtype=float)
dpr_bbt[dpr_bbt<0]= np.nan

dpr_bbh=np.array(dpr[scan]['CSF']['heightBB'], dtype=float)
dpr_bbh[dpr_bbh<0]= np.nan
dpr_bbw=np.array(dpr[scan]['CSF']['widthBB'], dtype=float)
dpr_bbw[dpr_bbw<0]= np.nan

lat_ppi =  50.730519999999999,
lon_ppi = 7.071663

#for iii in range(len(dpr_pp[1,1,:])):
#    plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp[:,:,iii]))
#    plt.xlim(2,12)
#    plt.ylim(49,53)
plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp_surf))
plt.colorbar()

plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
plt.grid()
plt.xlim(2,12)
plt.ylim(49,53)
#plt.show()

#position = np.where((dpr_lat<51.) & (dpr_lat>50.0) & (dpr_lon < 8) & (dpr_lon > 6)  )
position = np.where((dpr_lat<50.95) & (dpr_lat>50.70) & (dpr_lon < 7.5) & (dpr_lon > 6.5)  )


lat = dpr_lat[position]
lon = dpr_lon[position]
pp = dpr_pp[position]
dpr_time = dpr['NS']['ScanTime']

stunde = np.array(dpr_time['Hour'])[position[0]]
minute = np.array(dpr_time['Minute'])[position[0]]
sekunde = np.array(dpr_time['Second'])[position[0]]

jahr = np.array(dpr_time['Year'])[position[0]]
monat = np.array(dpr_time['Month'])[position[0]]
tag = np.array(dpr_time['DayOfMonth'])[position[0]]

zeit =  (str(jahr)+'.'+str(monat)+'.'+str(tag) + ' -- ' + str(stunde)+':'+str(minute)+':'+str(sekunde))

hdpr = 1000 * (np.arange(176,0,-1)*0.125) # Bei 88 500m und bei 176 ist es 250m

ff = 20

plt.subplot(2,2,2)

for jjj in range(len(pp[:,1])):
    plt.plot(pp[jjj,:], hdpr)

#plt.title('DPR '+ zeit, fontsize=ff)
plt.xlabel('Reflectivity in dBZ', fontsize=ff)
plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.ylim(0,6000)
plt.xlim(0,50)

plt.subplot(2,2,3)
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbh), vmin=0, vmax=4000)
plt.colorbar()

plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
plt.grid()
plt.xlim(2,12)
plt.ylim(49,53)
plt.subplot(2,2,4)
plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbw), vmin=0, vmax=4000)
plt.colorbar()

plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
plt.grid()
plt.xlim(2,12)
plt.ylim(49,53)

plt.show()


#plt.plot(pp[0,:], hdpr)
#plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)
#plt.xlabel('Reflectivity in dBZ')
#plt.grid()
#plt.show()
