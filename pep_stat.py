
import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace

import matplotlib.pyplot as plt

para = ['RR_MS', 'RR_NS', 'RR_HS', 'REF_MS', 'REF_NS', 'REF_HS']

# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]

pp = 0

TH = thresh[pp]

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


    #if pfad_gpm[ii][56:58] in sommer:
    R = h5py.File(pfad_gpm[ii], 'r')
    rrr = np.array(R['dat_rad'])
    ggg = np.array(R['dat_sat'])
    #xx=  np.array(R['x'])
    #yy = np.array(R['y'])
    rrr[rrr<=0]=np.nan
    ggg[ggg<=0]=np.nan

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
    #else:
    #    pass
    #print r_rad.shape

print r_rad

r_rad[r_rad<=TH]=np.nan
r_sat[r_sat<=TH]=np.nan


print r_rad.shape, r_sat.shape
print '------RADAR ------------ Satellite'
print 'Min: ',np.nanmin(r_rad), np.nanmin(r_sat)
print 'Mean: ', np.nanmean(r_rad), np.nanmean(r_sat)
print 'Max: ',np.nanmax(r_rad), np.nanmax(r_sat)

maske = ~np.isnan(r_sat) & ~np.isnan(r_rad)
slope, intercept, r_value, p_value, std_err = stats.linregress(r_sat[maske], r_rad[maske])
line = slope * r_sat +intercept

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(r_rad[maske], r_sat[maske])
line2 = slope2 * r_rad +intercept2

diffi = r_sat[maske] - r_rad[maske]
bias = np.nansum(diffi)/len(diffi)

biaslog = 10*np.log10(np.nansum(r_sat[maske])/np.nansum(r_rad[maske]))

ff = 20
cc = 1
fig = plt.figure(figsize=(12,12))
ax4 = fig.add_subplot(111, aspect='equal')
ax4.scatter(r_sat, r_rad, color='grey', alpha=0.6)

r_value_s, p_value_s = stats.spearmanr(r_sat[maske],r_rad[maske])

text = (#'f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
           '\n  Corr: ' + str(round(r_value,3)) + r'$\pm$ '+  str(round(std_err,3))+
        '\n  bias: '+ str(round(bias,3))+
        '\n  bias_log: '+ str(round(biaslog,3)))

ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=20,
                horizontalalignment='left', verticalalignment='top')

t1 = linspace(0,50,50)
plt.plot(t1,t1,'k--')
plt.plot(t1,t1*slope + intercept,label='RADOLAN', color='green', lw=1.5)
plt.plot(t1*slope2 + intercept2,t1,label='GPM DPR', color='blue', lw=1.5)
#plt.plot(t1,t1*slope + intercept , 'r-', lw=3 ,label='Regression')
#plt.plot(t1, t1*slope2 + intercept2, 'r-', lw=3 ,label='Regression')


plt.legend(loc='lower right', fontsize=15, scatterpoints= 1, numpoints=1,
           shadow=True, title='lin. Regression. (Ref.)')

#plt.xlim(0,50)
#plt.ylim(0,50)


plt.xlabel('GPM DPR Reflectivity (dBZ)',fontsize=ff)
plt.ylabel('RADOLAN Reflectivity (dBZ)',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid(color='r')
plt.title(str(para[pp]))
plt.tight_layout()



plt.show()


from satlib import validation_plot
validation_plot(r_sat, r_rad, TH)
plt.title(str(para[pp]))
plt.show()

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