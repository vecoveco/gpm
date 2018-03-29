import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr


import matplotlib.pyplot as plt
from pcc import get_my_cmap2

# Threshold ATBD GPM 2016
#thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]
thresh = [0.2, 0.5, 0.2, 15.0, 15.0, 15.0]



def corcorcor(A,B):
    mask = ~np.isnan(A) & ~np.isnan(B)
    diffi = A[mask]- B[mask]
    bias = np.nansum(diffi)/len(diffi)
    slope, intercept, r_value, p_value, std_err = stats.linregress(A[mask], B[mask])
    line = slope * A +intercept
    return r_value, std_err, bias


para = ['NS', 'HS','MS']
freq = ['13.6 GHz','35.5 GHz','35.5 GHz']
band = ['Ku','Ka','Ka']

pp = 0
dataS = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/SS_'+str(para[pp])+'.npy')

dataW = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/WS_'+str(para[pp])+'.npy')

#dataall = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/S_'+str(para[pp])+'.npy')

data = dataW
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

ff, ff2 = 20, 20



xbins = 10**np.arange(-1, 2.08, 0.07)
ybins = 10**np.arange(-1, 2.08, 0.07)

counts, _, _ = np.histogram2d(A[maske_p],B[maske_p], bins=(xbins, ybins))
'''
fig, ax = plt.subplots(figsize=(8,8))
#ax1 = fig.add_subplot(221, aspect='equal')
plt.pcolormesh(xbins, ybins, counts.T, cmap=get_my_cmap2(),vmin=0.1,vmax=600 )
cb = plt.colorbar()
cb.set_label('number of samples', fontsize=ff)

ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('DPR Rainrate in mm/h',fontsize=ff)
plt.ylabel('Radolan Rainrate in mm/h',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.title('Winter Corr: ' + str(round(corr,3)) + r'$\pm$'+  str(round(eror,3)), fontsize=20)

plt.grid()
plt.show()

#plt.savefig('/automount/ags/velibor/plot/IRS/NSwinter_RR.png')
#plt.close()
'''


fig =plt.figure(figsize=(8,8))
ax2 = fig.add_subplot(111, aspect='equal')
plt.hist2d(C[maske_z],D[maske_z], bins=60, cmap=get_my_cmap2(),vmin=0.1, vmax=400)
cbar = plt.colorbar(shrink=0.7)
cbar.set_label('number of samples', fontsize=ff2)

#cbar.set_ticks(labelsize=20)
cbar.ax.tick_params(labelsize=ff)

cx,cy = np.arange(0,80,1),np.arange(0,80,1)
plt.plot(cx,cy, color='black')
plt.title('Winter Corr: ' + str(round(corr2,3)) + r'$\pm$'+  str(round(eror2,3)), fontsize=ff2)
plt.xlim(15,70)
plt.ylim(15,70)

plt.xlabel('DPR Reflectivity in dBZ',fontsize=ff2)
plt.ylabel('Radolan Reflectivity in dBZ',fontsize=ff2)

plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

#plt.savefig('/automount/ags/velibor/plot/IRS/NSwinter_Ref.png')
#plt.close()
plt.show()
