"""

Vergleich der Radar Momente von  RADOLAN und deren ZEITLICH VERSCHIEBUNG

"""


import h5py
import numpy as np
import glob
import wradlib
import datetime as dt
from osgeo import osr
import matplotlib.pyplot as plt
from satlib import read_rado

t1, t2 = "20170307024500", "20170307025000"

x, y, z, bz = read_rado(t1, 'rx')
x1, y1, z1, bz1 = read_rado(t2, 'rx')

from pcc import get_my_cmap
import pcc
from scipy import stats, linspace


ff = 15
cc = 0.5
fig = plt.figure(figsize=(12,12))
ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------

pm1 = plt.pcolormesh(x, y, z, cmap=get_my_cmap(), vmin=0.01, vmax=50, zorder=2)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

pcc.plot_borders(ax1)

plt.title('RADOLAN Reflectivity:\n '+ t1+' UTC',fontsize=ff)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-400,450)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(222, aspect='equal')#------------------------------------

pm2 = plt.pcolormesh(x, y, z1, cmap=get_my_cmap(), vmin=0.01, vmax=50, zorder=2)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

pcc.plot_borders(ax2)

plt.title('RADOLAN Reflectivity:\n '+t2+' UTC',fontsize=ff)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-400,450)
plt.ylim(-4700, -3700)


ax3 = fig.add_subplot(223, aspect='equal')#------------------------------------

a, a1 = z.copy(), z1.copy()

b = a- a1

#a[a<0]=np.nan
#a1[a1<0]=np.nan
b[(a<0) & (a1<0)]=np.nan

print(np.nanmean(abs(b)))

pm3 = plt.pcolormesh(x, y,np.ma.masked_invalid(b), cmap='bwr', vmin=-15, vmax=15, zorder=2)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

pcc.plot_borders(ax3)

plt.title('RADOLAN Reflectivity Difference',fontsize=ff)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-400,450)
plt.ylim(-4700, -3700)


ax4 = fig.add_subplot(224, aspect='equal')#------------------------------------

rrr, ggg = z.copy(), z1.copy()
rrr[rrr<0]=np.nan
ggg[ggg<0]=np.nan

maske = ~np.isnan(ggg) & ~np.isnan(rrr)



slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
line = slope * ggg +intercept

from pcc import skill_score
SS = skill_score(ggg,rrr,th=0)

ax4.scatter(ggg, rrr, label='Reflectivity [dBZ]', color='grey', alpha=0.6)

r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])

text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
           '\nCorr: ' + str(round(r_value,3)) + r'$\pm$ '+  str(round(std_err,3))+
        '\nN: '+ str(int(SS['N']))+
        '\nHit: ' + str(SS['H'])+
        '\nMiss: ' + str(SS['M'])+
        '\nFalse: ' + str(SS['F'])+
        '\nCnegative: ' + str(SS['C'])+
        '\nHR: ' + str(round(SS['HR'],3))+
        '\nPOD: ' + str(round(SS['POD'],3))+
        '\nFAR: ' + str(round(SS['FAR'],3))+
        '\nBID: ' + str(round(SS['BID'],3))+
        '\nHSS: ' + str(round(SS['HSS'],3))+
        '\nBias: '+ str(round(SS['bias'],3))+
        '\nRMSE: '+ str(round(SS['RMSE'],3))+
        '\nCorrS:' +  str(round(r_value_s,3))
        )

ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top')

t1 = linspace(0,50,50)
plt.plot(t1,t1,'k-')
plt.plot(t1, t1*slope + intercept, 'r-', lw=3 ,label='Regression')
plt.plot(t1, t1*slope + (intercept+5), 'r-.', lw=1.5 ,label=r'Reg $\pm$ 5 mdBZ')
plt.plot(t1, t1*slope + (intercept-5), 'r-.', lw=1.5 )
#plt.plot(np.nanmean(ggg),np.nanmean(rrr), 'ob', lw = 4,label='Mean')
#plt.plot(np.nanmedian(ggg),np.nanmedian(rrr), 'vb', lw = 4,label='Median')

import matplotlib as mpl
mean = [ np.nanmean(ggg),np.nanmean(rrr)]
width = np.nanstd(ggg)
height = np.nanstd(rrr)
angle = 0
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                          angle=180+angle, color='blue', alpha=0.8,
                          fill=False, ls='--', label='Std')
#ax4.add_patch(ell)

plt.legend(loc='lower right', fontsize=10, scatterpoints= 1, numpoints=1, shadow=True)

plt.xlim(0,50)
plt.ylim(0,50)


plt.xlabel('RADOLAN(t_0) Reflectivity [dBZ]',fontsize=ff)
plt.ylabel('RADOLAN(t_1) Reflectivity [dBZ]',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid(color='r')
plt.show()




from satlib import cp_dist as cp

cp(z, z1)
#b[np.isnan(b)]= -9999

#plt.hist(b)
#plt.show()