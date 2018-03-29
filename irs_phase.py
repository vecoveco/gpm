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
#data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/WS_'+str(para[pp])+'.npy')
#data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/SS_'+str(para[pp])+'.npy')
data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/all'+str(para[pp])+'.npy')

#gx = data[0,:]
#gy = data[1,:]

#g_bbh = data[2,:]
#g_bbh[g_bbh<=0]=np.nan

#g_bbw = data[3,:]
#g_bbw[g_bbw<=0]=np.nan

#g_top = data[4,:]
g_type = data[5,:]
g_phase = data[6,:]

g_p2d = data[7,:]
g_p2d[g_p2d<=0]=np.nan

g_z2d = data[8,:]
g_z2d[g_z2d<=0]=np.nan

g_ry = data[9,:]
g_ry[g_ry<=0]=np.nan

g_rx = data[10,:]
g_rx[g_rx<=0]=np.nan


#Threshold bestimmen
#A,B = g_p2d.copy(), g_ry.copy()
#A[A<0.1]=np.nan; B[B<0.1]=np.nan
#C, D = g_z2d.copy(), g_rx.copy()
#C[C<15]=np.nan; D[D<15]=np.nan
a, b = g_p2d[g_phase == 2].copy(), g_ry[g_phase == 2].copy()
c, d = g_p2d[g_phase == 0].copy(), g_ry[g_phase == 0].copy()

maske_ab = ~np.isnan(a) & ~np.isnan(b)
maske_cd = ~np.isnan(c) & ~np.isnan(d)

th = 0.1
a[a<=th]=np.nan; b[b<=th]=np.nan; c[c<=th]=np.nan; d[d<=th]=np.nan

r_corr,r_eror,r_bias = corcorcor(a,b)
r_corr2,r_eror2,r_bias2 = corcorcor(c,d)


A, B = g_z2d[g_phase == 2].copy(), g_rx[g_phase == 2].copy()
C, D = g_z2d[g_phase == 0].copy(), g_rx[g_phase == 0].copy()

th = 15
A[A<=th]=np.nan; B[B<=th]=np.nan; C[C<=th]=np.nan; D[D<=th]=np.nan

maske_AB = ~np.isnan(A) & ~np.isnan(B)
maske_CD = ~np.isnan(C) & ~np.isnan(D)

z_corr,z_eror,z_bias = corcorcor(A,B)
z_corr2,z_eror2,z_bias2 = corcorcor(C,D)


ff, ff2 = 20, 20

'''
fig =plt.figure(figsize=(6,14))
ax1 = fig.add_subplot(2, 1, 1)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.scatter(a,b,color='blue',alpha=0.3, label='Liquid Corr: ' + str(round(r_corr,3)) + r'$\pm$'+  str(round(r_eror,3)))#+
            #'\nBias: '+ str(round(r_bias,3)))
ax1.scatter(c,d,color='black',alpha=0.3, label='Solid Corr: ' + str(round(r_corr2,3)) + r'$\pm$'+  str(round(r_eror2,3)))#+
            #'\nBias: '+ str(round(r_bias2,3)))
#plt.hist2d(g_p2d[maske_p],g_ry[maske_p],bins=99, cmap='jet')
#plt.colorbar()
plt.legend(loc= 'upper left', scatterpoints= 1, fontsize=ff)
plt.xlim(0,300)
plt.ylim(0,300)

plt.xlabel('DPR Rainrate in mm/h',fontsize=ff)
plt.ylabel('Radolan Rainrate in mm/h',fontsize=ff)
#plt.title(para[pp] + ' Rainrates: \n'
#          'RADOLAN and DPR ',fontsize=ff2)# + band[pp]+ ' [' +freq[pp]+ ']')
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid()


ax2 = fig.add_subplot(2, 1, 2)
plt.scatter(A,B,color='green',alpha=0.3, label='Liquid Corr: ' + str(round(z_corr,3)) + r'$\pm$'+ str(round(z_eror,3)))#+
            #'\nBias: '+ str(round(z_bias,3)))
plt.scatter(C,D,color='black',alpha=0.3, label='Solid Corr: ' + str(round(z_corr2,3)) + r'$\pm$'+ str(round(z_eror2,3)))#+
            #'\nBias: '+ str(round(z_bias2,3)))
#plt.hist2d(g_z2d[maske_z],g_rx[maske_z],bins=99, cmap='jet')
#plt.colorbar()
plt.legend(loc='upper left', scatterpoints= 1, fontsize=ff)
plt.xlim(15,70)
plt.ylim(15,70)

plt.xlabel('DPR Reflectivity in dbz',fontsize=ff)
plt.ylabel('Radolan Reflectivity in dbz',fontsize=ff)
#plt.title(para[pp] + ' Reflectivities: \n'
#          'RADOLAN and DPR ',fontsize=ff2)# + band[pp]+ ' [' +freq[pp]+ ']')
plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)


fig.subplots_adjust(bottom=0.05, top=0.93,hspace=0.3, left=0.15)
#plt.tight_layout()
plt.show()
'''



'''
a, b = g_z2d[g_phase == 2], g_rx[g_phase == 2]
c, d = g_z2d[g_phase == 0], g_rx[g_phase == 0]

th = 15
a[a<=th]=np.nan; b[b<=th]=np.nan; c[c<=th]=np.nan; d[d<=th]=np.nan

fig =plt.figure(figsize=(14,6))
maske = ~np.isnan(a) & ~np.isnan(b)
maske2 = ~np.isnan(c) & ~np.isnan(d)

ax1 = fig.add_subplot(1, 1, 1)
ax1.scatter(a,b,color='blue',alpha=0.3, label='liquid: '+ corcor(a,b))
ax1.scatter(c,d,color='black',alpha=0.3, label='solid: '+ corcor(c,d))
ax1.legend()
plt.xlabel('dpr reflectivity in dbz')
plt.ylabel('radolan reflectivity in dbz')
plt.grid()
plt.show()'''

from pcc import get_my_cmap2

xbins = 10**np.arange(-1, 2.08, 0.08)
ybins = 10**np.arange(-1, 2.08, 0.08)

counts, _, _ = np.histogram2d(c[maske_cd],d[maske_cd], bins=(xbins, ybins))
'''
fig, ax = plt.subplots(figsize=(8,8))
plt.pcolormesh(xbins, ybins, counts.T, cmap=get_my_cmap2(),vmin=0.1)
cb = plt.colorbar()
cb.set_label('number of samples', fontsize=ff)

ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('DPR Rainrate in mm/h',fontsize=ff)
plt.ylabel('Radolan Rainrate in mm/h',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.title('Solid Corr: ' + str(round(r_corr2,3)) + r'$\pm$'+  str(round(r_eror2,3)), fontsize=20)

plt.grid()
#plt.savefig('/automount/ags/velibor/plot/IRS/NSsolidWS_RR.png')
#plt.close()
plt.show()
'''

from pcc import get_my_cmap2

fig =plt.figure(figsize=(8,8))
ax2 = fig.add_subplot(111, aspect='equal')
plt.hist2d(C[maske_CD],D[maske_CD], bins=60, cmap=get_my_cmap2(),vmin=0.1)
cbar = plt.colorbar(shrink=0.7)
cbar.set_label('number of samples', fontsize=ff2)
cbar.ax.tick_params(labelsize=ff)

plt.title('Solid Phase Corr: ' + str(round(z_corr2,3)) + r'$\pm$'+  str(round(z_eror2,3)), fontsize=ff2)
plt.xlim(15,70)
plt.ylim(15,70)
cx,cy = np.arange(0,80,1),np.arange(0,80,1)
plt.plot(cx,cy, color='black')
plt.xlabel('DPR Reflectivity in dBZ',fontsize=ff2)
plt.ylabel('Radolan Reflectivity in dBZ',fontsize=ff2)

plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
#plt.savefig('/automount/ags/velibor/plot/IRS/NSsolidWS_Ref.png')
#plt.close()
plt.show()












############ LIQUID

from pcc import get_my_cmap2

xbins = 10**np.arange(-1, 2.08, 0.08)
ybins = 10**np.arange(-1, 2.08, 0.08)

counts, _, _ = np.histogram2d(a[maske_ab],b[maske_ab], bins=(xbins, ybins))
'''
fig, ax = plt.subplots(figsize=(8,8))
plt.pcolormesh(xbins, ybins, counts.T, cmap=get_my_cmap2(),vmin=0.1)
cb = plt.colorbar()
cb.set_label('number of samples', fontsize=ff)

ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('DPR Rainrate in mm/h',fontsize=ff)
plt.ylabel('Radolan Rainrate in mm/h',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)

plt.title('Liquid Corr: ' + str(round(r_corr,3)) + r'$\pm$'+  str(round(r_eror,3)), fontsize=20)

plt.grid()
#plt.savefig('/automount/ags/velibor/plot/IRS/NSliquidWS_RR.png')
#plt.close()
plt.show()
'''

from pcc import get_my_cmap2

fig =plt.figure(figsize=(8,8))
ax2 = fig.add_subplot(111, aspect='equal')
plt.hist2d(A[maske_AB],B[maske_AB], bins=60, cmap=get_my_cmap2(),vmin=0.1)
cbar = plt.colorbar(shrink=0.7)
cbar.set_label('number of samples', fontsize=ff2)
cbar.ax.tick_params(labelsize=ff)

plt.title('Liquid Phase Corr: ' + str(round(z_corr,3)) + r'$\pm$'+  str(round(z_eror,3)), fontsize=ff2)
plt.xlim(15,70)
plt.ylim(15,70)
cx,cy = np.arange(0,80,1),np.arange(0,80,1)
plt.plot(cx,cy, color='black')
plt.xlabel('DPR Reflectivity in dBZ',fontsize=ff2)
plt.ylabel('Radolan Reflectivity in dBZ',fontsize=ff2)

plt.grid()
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
#plt.savefig('/automount/ags/velibor/plot/IRS/NSliquidWS_Ref.png')
#plt.close()
plt.show()