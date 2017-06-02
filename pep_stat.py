
import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace

import matplotlib.pyplot as plt

pfad = ('/automount/ags/velibor/gpmdata/dumpdata/RR_NS*.hdf5')
pfad_gpm = sorted(glob.glob(pfad))


r_rad, r_sat = [], []

for ii in range(len(pfad_gpm)):
    print ii
    R = h5py.File(pfad_gpm[ii], 'r')
    rrr = np.array(R['dat_rad'])
    ggg = np.array(R['dat_sat'])
    #x=  np.array(R['x'])
    #y = np.array(R['y'])
    rrr = rrr.reshape(rrr.shape[0]*rrr.shape[1])
    ggg = ggg.reshape(ggg.shape[0]*ggg.shape[1])
    #print rrr.shape


    r_rad = np.concatenate((r_rad, rrr), axis=0)
    r_sat = np.concatenate((r_sat, ggg), axis=0)

    #print r_rad.shape


TH = 0.5
r_rad[r_rad<=TH]=np.nan
r_sat[r_sat<=TH]=np.nan


print r_rad.shape, r_sat.shape
print np.nanmin(r_rad), np.nanmin(r_sat)
print np.nanmean(r_rad), np.nanmean(r_sat)
print np.nanmax(r_rad), np.nanmax(r_sat)

maske = ~np.isnan(r_sat) & ~np.isnan(r_rad)
slope, intercept, r_value, p_value, std_err = stats.linregress(r_sat[maske], r_rad[maske])
line = slope * r_sat +intercept

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(r_rad[maske], r_sat[maske])
line2 = slope2 * r_rad +intercept2

diffi = r_sat[maske]-r_rad[maske]
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

plt.tight_layout()



plt.show()


from satlib import validation_plot
validation_plot(r_sat, r_rad, TH)