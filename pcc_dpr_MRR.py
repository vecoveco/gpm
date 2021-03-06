"""

Reading MRR Data

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO

# MRR
pfad_mrr = '/automount/mrr/mrr2/2014/2014-10/2014-10-07/AveData_mrr2_20141007023828.ave.gz'
#pfad_mrr = '/automount/mrr/mrr2/2016/2016-01/2016-01-07/AveData_mrr2_20160107124312.ave.gz'


# DPR
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20141007-S015721-E032951.003445.V04A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20160107-S120629-E133900.010562.V04A.HDF5'

#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170306-S151005-E164237.017160.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170223-S182621-E195854.016991.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170222-S113453-E130727.016971.V05A.HDF5'
#dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170401-S003803-E021037.017555.V05A.HDF5'
dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170603-S235100-E012332.018550.V05A.HDF5'

df = pd.read_csv(pfad_mrr, compression='gzip', header=1, delim_whitespace=True,index_col=False)
df = df.set_index(u'H')


h = np.arange(150,4800,150)


plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)
plt.plot(df.loc['PIA'].values,h,label='PIA in dB', color='blue', linestyle='-.', lw=2)
plt.plot(df.loc['z'].values,h,label='att. Ref in dBZ', color='blue', linestyle='--', lw=2)
plt.plot(df.loc['TF'].values,h,label='TF', color='grey')
plt.plot(df.loc['RR'].values,h,label='RR', linestyle='-', color='black',lw=2)
plt.plot(df.loc['LWC'].values,h,label='Liquid Water Content')
plt.plot(df.loc['W'].values,h,label='Fallgeschwindigkeit')
plt.grid()
plt.legend(loc='lower right')
plt.xlabel('Reflectivity in dBZ')
plt.ylabel('Hight in m')
plt.title('MRR - ' + pfad_mrr[44:44+28])
plt.ylim(0,6000)
plt.xlim(0,50)
plt.show()




import satlib as sl
import h5py

scan  = 'NS' #or MS

dpr = h5py.File(dpr_pfad, 'r')
dpr_lat=np.array(dpr[scan]['Latitude'])
dpr_lon=np.array(dpr[scan]['Longitude'])
dpr_pp=np.array(dpr[scan]['SLV']['zFactorCorrected'])
dpr_pp[dpr_pp<0]= np.nan

lat_ppi =  50.730519999999999,
lon_ppi = 7.071663

for iii in range(len(dpr_pp[1,1,:])):
    plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp[:,:,iii]))
    plt.xlim(2,12)
    plt.ylim(49,53)

plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
plt.grid()
plt.show()

position = np.where((dpr_lat<51.) & (dpr_lat>50.0) & (dpr_lon < 8) & (dpr_lon > 6)  )
#position = np.where((dpr_lat<50.95) & (dpr_lat>50.70) & (dpr_lon < 7.5) & (dpr_lon > 6.5)  )


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

print zeit

hdpr = 1000 * (np.arange(176,0,-1)*0.125) # Bei 88 500m und bei 176 ist es 250m

ff = 20

plt.figure(figsize=(12,12))
plt.subplot(1,2,1)
plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)

plt.plot(df.loc['PIA'].values,h,label='PIA in dB', color='blue', linestyle='-.', lw=2)
plt.plot(df.loc['z'].values,h,label='att. Ref in dBZ', color='blue', linestyle='--', lw=2)
#plt.plot(df.loc['TF'].values,h,label='TF', color='grey')
plt.plot(df.loc['RR'].values,h,label='RR', linestyle='-', color='black',lw=2)
#plt.plot(df.loc['LWC'].values,h,label='Liquid Water Content')
#plt.plot(df.loc['W'].values,h,label='Fallgeschwindigkeit')
plt.grid()
plt.legend(loc='lower right', fontsize=ff)
plt.ylabel('Hight in m', fontsize=ff)
plt.title('MRR - ' + pfad_mrr[44:44+28], fontsize=ff)
plt.ylim(0,6000)
plt.xlim(0,50)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.subplot(1,2,2)
for jjj in range(len(pp[:,1])):
    plt.plot(pp[jjj,:], hdpr)

plt.title('DPR '+ zeit, fontsize=ff)
plt.xlabel('Reflectivity in dBZ', fontsize=ff)
plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.ylim(0,6000)
plt.xlim(0,50)
plt.show()


#plt.plot(pp[0,:], hdpr)
#plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)
#plt.xlabel('Reflectivity in dBZ')
#plt.grid()
#plt.show()
