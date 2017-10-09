import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr


import matplotlib.pyplot as plt



para = ['NS', 'HS','MS']



# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]


'''
for ip in range(3):

    pp = ip

    #pfad = ('/automount/ags/velibor/gpmdata/dumpdata/' + para[pp] + '/' + para[pp] + '*.hdf5')
    pfad = ('/automount/ags/velibor/gpmdata/dumpdata/gpm_dpr_rado_all2/dprrado_'+ para[pp]+'/dprrado'+'*.hdf5')

    pfad_gpm = sorted(glob.glob(pfad))

    print 'Es sind ',len(pfad_gpm), 'Dateien vorhanden!'

    gx, gy = [], []
    g_bbh, g_bbw, g_top = [], [], []
    g_type,g_phase = [], []
    g_p2d, g_z2d = [], []
    g_rx, g_ry = [], []


    for ii in range(len(pfad_gpm)):

        print ii,' von ', len(pfad_gpm)

        #if (pfad_gpm[ii][84:86] == '06') or (pfad_gpm[ii][84:86] == '07') or (pfad_gpm[ii][84:86] == '08'):
        if (pfad_gpm[ii][84:86] == '12') or (pfad_gpm[ii][84:86] == '01') or (pfad_gpm[ii][84:86] == '02'):
            print pfad_gpm[ii][84:86]

            R = h5py.File(pfad_gpm[ii], 'r')
            s1, s2 = 0, 48
            x = np.array(R['x'])[:,s1:s2]
            y = np.array(R['y'])[:,s1:s2]
            sat_bbh = np.array(R['sat_bbh'])[:,s1:s2]
            sat_bbw = np.array(R['sat_bbw'])[:,s1:s2]
            sat_phase = np.array(R['sat_phase'])[:,s1:s2]
            sat_type = np.array(R['sat_type'])[:,s1:s2]
            sat_p2d = np.array(R['sat_p2d'])[:,s1:s2]
            sat_z2d = np.array(R['sat_z2d'])[:,s1:s2]
            sat_top = np.array(R['sat_top'])[:,s1:s2]
            radolan_z2d = np.array(R['radolan_z2d'])[:,s1:s2]
            radolan_p2d = np.array(R['radolan_p2d'])[:,s1:s2]


            x = x.reshape(x.shape[0] * x.shape[1])
            y = y.reshape(y.shape[0] * y.shape[1])
            bbh = sat_bbh.reshape(sat_bbh.shape[0] * sat_bbh.shape[1])
            bbw = sat_bbw.reshape(sat_bbw.shape[0] * sat_bbw.shape[1])
            phase = sat_phase.reshape(sat_phase.shape[0] * sat_phase.shape[1])
            type = sat_type.reshape(sat_type.shape[0] * sat_type.shape[1])
            p2d = sat_p2d.reshape(sat_p2d.shape[0] * sat_p2d.shape[1])
            z2d = sat_z2d.reshape(sat_z2d.shape[0] * sat_z2d.shape[1])
            top = sat_top.reshape(sat_top.shape[0] * sat_top.shape[1])
            ry = radolan_p2d.reshape(radolan_p2d.shape[0] * radolan_p2d.shape[1])
            rx = radolan_z2d.reshape(radolan_z2d.shape[0] * radolan_z2d.shape[1])

            gx = np.concatenate((gx, x), axis=0)
            gy = np.concatenate((gy, y), axis=0)
            g_bbh = np.concatenate((g_bbh, bbh), axis=0)
            g_bbw = np.concatenate((g_bbw, bbw), axis=0)
            g_phase = np.concatenate((g_phase, phase), axis=0)
            g_type = np.concatenate((g_type, type), axis=0)
            g_p2d = np.concatenate((g_p2d, p2d), axis=0)
            g_z2d = np.concatenate((g_z2d, z2d), axis=0)
            g_top = np.concatenate((g_top, top), axis=0)
            g_ry = np.concatenate((g_ry, ry), axis=0)
            g_rx = np.concatenate((g_rx, rx), axis=0)

            g_ry[g_ry < 0.01] = np.nan
            g_rx[g_rx < 0.01] = np.nan

        else:
            print 'nicht in season ', pfad_gpm[ii][84:86]

    np.save('/automount/ags/velibor/gpmdata/dumpdata/npy/WS_'+str(para[pp]),
            (gx,gy,g_bbh, g_bbw, g_top, g_type, g_phase, g_p2d,
             g_z2d, g_ry, g_rx))

'''


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
pp = 2
data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/WS_'+str(para[pp])+'.npy')
#data = np.load('/automount/ags/velibor/gpmdata/dumpdata/npy/WS_'+str(para[pp])+'.npy')

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

ff, ff2 = 15, 20

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


