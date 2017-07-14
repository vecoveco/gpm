import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace

import matplotlib.pyplot as plt

para = ['RR_MS', 'RR_NS', 'RR_HS', 'REF_MS', 'REF_NS', 'REF_HS']

# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]

for jjj in range(3):
    pp = jjj

    #TH = thresh[pp]

    #pfad = ('/automount/ags/velibor/gpmdata/dumpdata/RR_MS/RR_MS*.hdf5')
    pfad = ('/automount/ags/velibor/gpmdata/dumpdata/' + para[pp] + '/' + para[pp] + '*.hdf5')

    pfad_gpm = sorted(glob.glob(pfad))

    print 'Es sind ',len(pfad_gpm), 'Dateien vorhanden!'

    sommer, winter, = ['07','08','06'], ['12','01','02']
    r_rad, r_sat = [], []
    rmax_rad, rmax_sat = [], []
    rmin_rad, rmin_sat = [], []
    rm_rad, rm_sat = [], []

    #len(pfad_gpm)

    for ii in range(len(pfad_gpm)):
        #print ii

        # by REF nicht 56:58 sonder 58:60
        if pfad_gpm[ii][56:58] in sommer:
            #print pfad_gpm[ii][56:58]
            R = h5py.File(pfad_gpm[ii], 'r')
            rrr = np.array(R['dat_rad'])
            ggg = np.array(R['dat_sat'])
            #xx=  np.array(R['x'])
            #yy = np.array(R['y'])
            rrr[rrr<=0] = np.nan
            ggg[ggg<=0] = np.nan

            rrr[rrr>1000000000] = np.nan
            ggg[ggg>1000000000] = np.nan

            # TODO: Chek einbauen ob der Uberflug relevant ist

            rrr = rrr.reshape(rrr.shape[0]*rrr.shape[1])
            ggg = ggg.reshape(ggg.shape[0]*ggg.shape[1])
            #print rrr.shape


            r_rad = np.concatenate((r_rad, rrr), axis=0)
            r_sat = np.concatenate((r_sat, ggg), axis=0)

            rmax_rad.append(np.nanmax(rrr))
            rmax_sat.append(np.nanmax(ggg))

            rmin_rad.append(np.nanmin(rrr))
            rmin_sat.append(np.nanmin(ggg))

            rm_rad.append(np.nanmean(rrr))
            rm_sat.append(np.nanmean(ggg))
        else:
            pass
        #print r_rad.shape
    print np.nanmin(r_rad), np.nanmin(r_sat)
    minirad, minisat = np.nanmin(r_rad), np.nanmin(r_sat)
    print minirad
    #print minirad, minisat

    maxrad, maxsat = np.nanmax(r_rad), np.nanmax(r_sat)
    #print maxrad, maxsat

    TH = np.nanmax(np.array([minirad, minisat]))
    #print TH
    #TH = 10

    r_rad[r_rad<=TH]=np.nan
    r_sat[r_sat<=TH]=np.nan


    #print r_rad.shape, r_sat.shape
    #print '------RADAR ------------ Satellite'
    #print 'Min: ',np.nanmin(r_rad), np.nanmin(r_sat)
    #print 'Mean: ', np.nanmean(r_rad), np.nanmean(r_sat)
    #print 'Max: ',np.nanmax(r_rad), np.nanmax(r_sat)


    from satlib import validation_plot, validation_plot_log
    if pp <= 2:
        print '-----------------------> RR'
        validation_plot_log(r_sat, r_rad, TH)
        plt.title(str(para[pp])+ '- TH: ' + str(TH))
        #plt.show()
        plt.savefig('/automount/ags/velibor/plot/validation/run/validation_dynTH_JJA'+ para[pp] +'.png' )
        plt.close()
    if pp>2:
        print '----------------------> REF'
        validation_plot(r_sat, r_rad, TH)
        plt.title(str(para[pp])+ '- TH: ' + str(TH))
        #plt.show()
        plt.savefig('/automount/ags/velibor/plot/validation/run/validation_dynTH_JJA'+ para[pp] +'.png' )
        plt.close()

'''
rm_rad = np.array(rm_rad)
rm_sat = np.array(rm_sat)
rmin_rad = np.array(rmin_rad)
rmin_sat = np.array(rmin_sat)
rmax_rad = np.array(rmax_rad)
rmax_sat = np.array(rmax_sat)

rm_rad[rm_rad<=0]=np.nan
rm_sat[rm_sat<=0]=np.nan
rmin_rad[rmin_rad<=0]=np.nan
rmin_sat[rmin_sat<=0]=np.nan
rmax_rad[rmax_rad<=0]=np.nan
rmax_sat[rmax_sat<=0]=np.nan

plt.subplot(3,1,1)
plt.plot(rm_rad)
plt.plot(rm_sat)

plt.subplot(3,1,2)
plt.plot(rmin_rad)
plt.plot(rmin_sat)

plt.subplot(3,1,3)
plt.plot(rmax_rad)
plt.plot(rmax_sat)

plt.title(str(para[pp]))

plt.show()
'''

maske = ~np.isnan(r_sat) & ~np.isnan(r_rad)

plt.boxplot([r_sat[maske], r_rad[maske]])
plt.show()