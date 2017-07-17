"""

Reading MRR Data

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO
import satlib as sl
import h5py
import glob

dates = ['20150128','20150330', '20150404','20151208','20151216','20160107','20160601','20160612','20161019','20161024','20161109','20161222',
         '20170113','20141007','20140708', '20140729']#, '20151015','20160209', '20160915','20161121', '20141008'

#dates =['20141007']

for ii in dates:
    zt=ii
    print zt
    pfad = ('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.'+zt+'*.HDF5')
    print pfad
    dpr_pfad = sorted(glob.glob(pfad))[0]
    print dpr_pfad

    # DPR

    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20141007-S015721-E032951.003445.V04A.HDF5'

    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20160107-S120629-E133900.010562.V04A.HDF5'
    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170306-S151005-E164237.017160.V05A.HDF5'
    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170223-S182621-E195854.016991.V05A.HDF5'
    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170222-S113453-E130727.016971.V05A.HDF5'
    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170401-S003803-E021037.017555.V05A.HDF5'
    #dpr_pfad = '/automount/ags/velibor/gpmdata/dpr_brandon_BB/2A.GPM.DPR.V7-20170308.20170603-S235100-E012332.018550.V05A.HDF5'







    scan  = 'NS' #or MS

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

    lat_ppi =  50.730519999999999,
    lon_ppi = 7.071663



    #position = np.where((dpr_lat<51.) & (dpr_lat>50.0) & (dpr_lon < 8) & (dpr_lon > 6)  )
    l1, l2 = 52.9, 49.01
    k1, k2 = 10.8, 3.5
    #l1, l2 = 50.90, 50.50
    #k1, k2 = 7.8, 6.5
    position = np.where((dpr_lat<l1) & (dpr_lat>l2) & (dpr_lon < k1) & (dpr_lon > k2)  )


    lat = dpr_lat[position]
    lon = dpr_lon[position]
    pp = dpr_pp[position]
    dpr_time = dpr['NS']['ScanTime']

    stunde = np.array(dpr_time['Hour'])[position[0]][0]
    minute = np.array(dpr_time['Minute'])[position[0]][0]
    sekunde = np.array(dpr_time['Second'])[position[0]][0]

    jahr = np.array(dpr_time['Year'])[position[0]][0]
    monat = np.array(dpr_time['Month'])[position[0]][0]
    tag = np.array(dpr_time['DayOfMonth'])[position[0]][0]

    zeit =  (str(jahr)+'.'+str(monat)+'.'+str(tag) + ' -- ' + str(stunde)+':'+str(minute)+':'+str(sekunde))

    print jahr, monat, tag, stunde, minute

    h = np.arange(150,4800,150)
    hdpr = 1000 * (np.arange(176,0,-1)*0.125)
    hhh = np.array(len(pp[:,0])*list(hdpr))
    ppp = pp.reshape(pp.shape[0]*pp.shape[1])

    maske = ~np.isnan(hhh) & ~np.isnan(ppp)

    fig = plt.figure(figsize=(12,12))
    zzz = str(jahr)+'-'+str(monat)+'-'+str(tag)+'--'+str(stunde)+':'+str(minute)+' UTC'
    fig.suptitle(zzz + ' UTC')
    plt.subplot(2,2,1)
    plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_pp_surf))
    plt.colorbar()

    ax1 = plt.scatter(lon_ppi, lat_ppi, c=50 ,s=50, color='red')

    plt.scatter(k1,l1, c=50 ,s=50, color='red')
    plt.scatter(k2,l1, c=50 ,s=50, color='red')
    plt.scatter(k1,l2, c=50 ,s=50, color='red')
    plt.scatter(k2,l2, c=50 ,s=50, color='red')
    print lon_ppi, lat_ppi

    plt.grid()
    plt.xlim(2,12)
    plt.ylim(49,53)
    plt.subplot(2,2,2)

    #for jjj in range(len(pp[:,1])):
    #    plt.plot(pp[jjj,:], hdpr, color='black')
    plt.hist2d(ppp[maske],hhh[maske], bins=30)
    plt.colorbar()


    plt.title('DPR Ref. in Box')
    plt.xlabel('Reflectivity in dBZ')
    plt.grid()
    plt.xticks()
    plt.yticks()

    plt.ylim(0,6000)
    plt.xlim(0,50)

    plt.subplot(2,2,3)
    plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbh), vmin=0, vmax=4000)
    plt.colorbar()

    plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
    plt.grid()
    plt.title('BB Hight')
    plt.xlim(2,12)
    plt.ylim(49,53)


    plt.subplot(2,2,4)
    plt.pcolormesh(dpr_lon, dpr_lat,np.ma.masked_invalid(dpr_bbw), vmin=0, vmax=1000)
    plt.colorbar()

    plt.scatter(lon_ppi, lat_ppi, c=100 ,s=100, color='red')
    plt.grid()
    plt.title('BB Width')
    plt.xlim(2,12)
    plt.ylim(49,53)

    #plt.show()
    plt.savefig('/automount/ags/velibor/plot/BB/'+'HS2dprbb_'+str(zzz)+'.png' )
    plt.close()

    #plt.plot(pp[0,:], hdpr)
    #plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)
    #plt.xlabel('Reflectivity in dBZ')
    #plt.grid()
    #plt.show()

