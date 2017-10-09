import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr


import matplotlib.pyplot as plt


def corcor(A,B):
    mask = ~np.isnan(A) & ~np.isnan(B)
    corr = np.corrcoef(A[mask],B[mask])[0,1]
    corr = round(corr,3)
    return str(corr)

def corcorcor(A,B):
    mask = ~np.isnan(A) & ~np.isnan(B)
    diffi = A[mask]-B[mask]
    bias = np.nansum(diffi)/len(diffi)
    slope, intercept, r_value, p_value, std_err = stats.linregress(A[mask], B[mask])
    line = slope * A +intercept
    return r_value, std_err, bias


para = ['NS', 'HS','MS']
freq = ['13.6 GHz','35.5 GHz','35.5 GHz']
band = ['Ku','Ka','Ka']
pp = 0
data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/all'+str(para[pp])+'.npy')

#gx = data[0,:]
#gy = data[1,:]

#g_bbh = data[2,:]
#g_bbh[g_bbh<=0]=np.nan

#g_bbw = data[3,:]
#g_bbw[g_bbw<=0]=np.nan

#g_top = data[4,:]
#g_type = data[5,:]
#g_phase = data[6,:]

g_p2d = data[7,:]
g_p2d[g_p2d<=0]=np.nan

g_z2d = data[8,:]
g_z2d[g_z2d<=0]=np.nan

g_ry = data[9,:]
g_ry[g_ry<=0]=np.nan

g_rx = data[10,:]
g_rx[g_rx<=0]=np.nan


#Threshold bestimmen
A,B = g_p2d.copy(), g_ry.copy()
A[A<0.1]=np.nan; B[B<0.1]=np.nan
C, D = g_z2d.copy(), g_rx.copy()
C[C<15]=np.nan; D[D<15]=np.nan


maske_p = ~np.isnan(A) & ~np.isnan(B)
maske_z = ~np.isnan(C) & ~np.isnan(D)

corr,eror,bias = corcorcor(A,B)
corr2,eror2,bias2 = corcorcor(C,D)

ff, ff2 =15, 20

fig =plt.figure(figsize=(6,14))
ax1 = fig.add_subplot(2, 1, 1)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.scatter(A,B,color='blue',alpha=0.3, label='Corr: ' + str(round(corr,3)) + r'$\pm$'+  str(round(eror,3))+
            '\nBias: '+ str(round(bias,3)))
#plt.hist2d(g_p2d[maske_p],g_ry[maske_p],bins=99, cmap='jet')
#plt.colorbar()
plt.legend(loc= 'upper left', scatterpoints= 1, fontsize=ff)
plt.xlim(0,300)
plt.ylim(0,300)

plt.xlabel('DPR Rainrate in mm/h',fontsize=ff)
plt.ylabel('Radolan Rainrate in mm/h',fontsize=ff)
plt.title(para[pp] + ' Rainrates: \n'
          'RADOLAN and DPR ',fontsize=ff2)# + band[pp]+ ' [' +freq[pp]+ ']')
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid()


ax2 = fig.add_subplot(2, 1, 2)
plt.scatter(C,D,color='green',alpha=0.3, label='Corr: ' + str(round(corr2,3)) + r'$\pm$'+ str(round(eror2,3))+
            '\nBias: '+ str(round(bias2,3)))
#plt.hist2d(g_z2d[maske_z],g_rx[maske_z],bins=99, cmap='jet')
#plt.colorbar()
plt.legend(loc='upper left', scatterpoints= 1, fontsize=ff)
plt.xlim(15,70)
plt.ylim(15,70)

plt.xlabel('DPR Reflectivity in dbz',fontsize=ff)
plt.ylabel('Radolan Reflectivity in dbz',fontsize=ff)
plt.title(para[pp] + ' Reflectivities: \n'
          'RADOLAN and DPR ',fontsize=ff2)# + band[pp]+ ' [' +freq[pp]+ ']')
plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)


fig.subplots_adjust(bottom=0.05, top=0.95,hspace=0.3, left=0.15)
#plt.tight_layout()
plt.show()


'''
plt.figure(figsize=(14,6))
plt.subplot(1,3,1)
plt.scatter(ns_z_dpr[ns_p_rad == 2],ns_z_rad[ns_p_rad == 2],
            label='liquid: '+ corcor(ns_z_dpr[ns_p_rad == 2],ns_z_rad[ns_p_rad == 2]))
plt.scatter(ns_z_dpr[ns_p_rad == 0],ns_z_rad[ns_p_rad == 0],
            label='solid: '+ corcor(ns_z_dpr[ns_p_rad == 0],ns_z_rad[ns_p_rad == 0]),color='black')
plt.legend()
plt.xlabel('dpr reflectivity in dbz')
plt.ylabel('radolan reflectivity in dbz')
plt.title('NS')
plt.grid()'''