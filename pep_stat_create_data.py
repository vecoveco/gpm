import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr


import matplotlib.pyplot as plt

#para = ['RR_MS', 'RR_NS', 'RR_HS', 'REF_MS', 'REF_NS', 'REF_HS', 'SKILL_NS', 'SKILL_HS', 'SKILL_MS']

para = ['NS', 'HS','MS']


#################################################################DEM
dem=np.load("/automount/ftp/velibor/dem.npy")
dem_x, dem_y,dem_v = dem[0,:,:],dem[1,:,:],dem[2,:,:]

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

dx, dy = wradlib.georef.reproject(dem_x, dem_y, projection_target=proj_stereo , projection_source=proj_wgs)

# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]

for ip in range(1):

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


    #len(pfad_gpm)

    for ii in range(len(pfad_gpm)):

        print ii,' von ', len(pfad_gpm)

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
        print np.nanmax(sat_phase)

        ############################################### DEM
        #gk3 = wradlib.georef.epsg_to_osr(31467)
        #grid_gpm_xy = np.vstack((x.ravel(), y.ravel())).transpose()
        #dxy = np.vstack((dx.ravel(), dy.ravel())).transpose()

        #dr = wrl.ipol.interpolate(dxy, grid_gpm_xy, dem_v.reshape(dem_v.shape[0]*dem_v.shape[1],1), wrl.ipol.Idw, nnearest=4)

        #dr = np.ma.masked_invalid(dr)

        #dem_h = dr.reshape(x.shape)

        #h = dem_h.reshape(dem_h.shape[0] * dem_h.shape[1])

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


    #np.save('/automount/ags/velibor/gpmdata/dumpdata/npy/2all'+str(para[pp]),
    #        (gx,gy,g_bbh, g_bbw, g_top, g_type, g_phase, g_p2d,
    #         g_z2d, g_ry, g_rx))

ggg = g_z2d.copy()
rrr = g_rx.copy()

ggg[ggg<15]=np.nan
rrr[rrr<15]=np.nan

ratio1 = ggg/rrr
#ratio2 = ggg-rrr
ratio2 = ggg/np.nanmean(rrr)
# Ratioof data

A, B = ratio1,ratio2
#A[A<15]=np.nan
#B[B<15]=np.nan
from satlib import corcor
maskr = ~np.isnan(A) & ~np.isnan(B)
plt.hist2d(A[maskr],B[maskr],bins=100,vmin=0.01)
plt.title(corcor(A[maskr],B[maskr]))
plt.xlabel('DPR/RADOLAN')
plt.ylabel('DPR - RADOLAN')
#plt.xlim(0,2)
plt.colorbar()
plt.show()
'''
ggg[g_phase != 2]=np.nan # liqud
rrr[g_phase != 2]=np.nan # liqud

#[g_phase == 0] # soild
AA, BB = A[maskr], B[maskr]
CC = g_type[maskr]

AAA, BBB, CCC = AA[(g_type[maskr] >=1) & (g_type[maskr] < 2.)], BB[(g_type[maskr] >=1) & (g_type[maskr] < 2.)], CC[(g_type[maskr] >=1) & (g_type[maskr] < 2.)]
DDD = [(g_type[maskr] >=2) & (g_type[maskr] < 3.)]
plt.scatter(AAA,
            BBB,
            c=DDD, alpha=0.2)
plt.colorbar()
plt.show()
'''


from satlib import corcor

plt.figure(1)
plt.scatter(g_z2d, g_p2d, color='green', alpha=0.5, label='dpr: '+corcor(g_z2d, g_p2d))
plt.scatter(g_rx, g_ry, color='blue', alpha=0.5, label='radolan: '+corcor(g_rx, g_ry))
plt.grid()
plt.legend()
plt.xlabel('Ref. in dbz')
plt.ylabel('Rain in mm/h')
plt.show()

plt.figure(2)
plt.subplot(1,2,1)
A, B = g_top, g_rx
maske = ~np.isnan(A) & ~np.isnan(B)
plt.hist2d(A[maske], B[maske], bins=99)
plt.colorbar()
plt.title(str(np.corrcoef(A[maske],B[maske])[0,1]))
plt.xlabel('dpr storm top in m')
plt.ylabel('radolan Ref. in dbz')
plt.grid()
plt.subplot(1,2,2)
A, B = g_top, g_z2d
maske = ~np.isnan(A) & ~np.isnan(B)
plt.hist2d(A[maske], B[maske], bins=99)
plt.colorbar()
plt.title(str(np.corrcoef(A[maske],B[maske])[0,1]))
plt.xlabel('dpr storm top in m')
plt.ylabel('dpr Ref. in dbz')
plt.grid()
plt.show()

plt.figure(3)
g_type = g_type/10000000
plt.scatter(g_z2d[(g_type >=1) & (g_type < 2.)], g_rx[(g_type >=1) & (g_type < 2.)],
            alpha=0.3, color='blue', label='stratiform dpr: '+corcor(g_z2d[(g_type >=1) & (g_type < 2.)], g_rx[(g_type >=1) & (g_type < 2.)]))

plt.scatter(g_z2d[(g_type >=2) & (g_type < 3.)], g_rx[(g_type >=2) & (g_type < 3.)], alpha=0.3,
            color='green', label='convektiv dpr: ' + corcor(g_z2d[(g_type >=2) & (g_type < 3.)], g_rx[(g_type >=2) & (g_type < 3.)]))

plt.scatter(g_z2d[(g_type >=3)], g_rx[(g_type >=3)], alpha=0.3,
            color='red', label='other dpr: ' + corcor(g_z2d[(g_type >=3)], g_rx[(g_type >=3)]))

plt.xlabel('DPR')
plt.ylabel('radolan')
plt.grid()
plt.legend()
plt.show()

plt.subplot(2,2,1)
plt.scatter(g_z2d[(g_type >=2) & (g_type < 3.) & (g_phase == 2)], g_rx[(g_type >=2) & (g_type < 3.)& (g_phase == 2)], alpha=0.3,
            color='green', label='convektiv dpr: ' + corcor(g_z2d[(g_type >=2) & (g_type < 3.)& (g_phase == 2)], g_rx[(g_type >=2) & (g_type < 3.)& (g_phase == 2)]))
plt.title('Convectiv and liquid')
plt.legend()
plt.grid()
plt.subplot(2,2,2)
plt.scatter(g_z2d[(g_type >=2) & (g_type < 3.) & (g_phase == 0)], g_rx[(g_type >=2) & (g_type < 3.)& (g_phase == 0)], alpha=0.3,
            color='green', label='convektiv dpr: ' + corcor(g_z2d[(g_type >=2) & (g_type < 3.)& (g_phase == 0)], g_rx[(g_type >=2) & (g_type < 3.)& (g_phase == 0)]))
plt.title('Convectiv and solid')
plt.legend()
plt.grid()
plt.subplot(2,2,3)
plt.scatter(g_z2d[(g_type >=1) & (g_type < 2.)& (g_phase == 2)], g_rx[(g_type >=1) & (g_type < 2.)& (g_phase == 2)],
            alpha=0.3, color='blue', label='stratiform dpr: '+corcor(g_z2d[(g_type >=1) & (g_type < 2.)& (g_phase == 2)], g_rx[(g_type >=1) & (g_type < 2.)& (g_phase == 2)]))
plt.title('stratiform and liquid')
plt.legend()
plt.grid()
plt.subplot(2,2,4)
plt.scatter(g_z2d[(g_type >=1) & (g_type < 2.)& (g_phase == 0)], g_rx[(g_type >=1) & (g_type < 2.)& (g_phase == 0)],
            alpha=0.3, color='blue', label='stratiform dpr: '+corcor(g_z2d[(g_type >=1) & (g_type < 2.)& (g_phase == 0)], g_rx[(g_type >=1) & (g_type < 2.)& (g_phase == 0)]))
plt.title('stratiform  and solid')
plt.legend()
plt.grid()
plt.show()

plt.figure(4)
plt.subplot(2,2,1)

plt.scatter(g_z2d[g_phase == 1], g_rx[g_phase == 1],
            alpha=0.2, color='blue', label='mixed')
plt.scatter(g_z2d[g_phase == 2], g_rx[g_phase == 2],
            alpha=0.2, color='red', label='liquid')
plt.scatter(g_z2d[g_phase == 0], g_rx[g_phase== 0],
            alpha=0.5, color='black', label='solid')
plt.xlabel('dpr')
plt.ylabel('radolan')

plt.legend()
plt.grid()

plt.subplot(2,2,2)
maske1 = ~np.isnan(g_z2d[g_phase == 0]) & ~np.isnan(g_rx[g_phase == 0])
plt.hist2d(g_z2d[g_phase  == 0][maske1], g_rx[g_phase == 0][maske1],bins=99)
plt.colorbar()
plt.title('solid - corr: '+ corcor(g_z2d[g_phase  == 0][maske1], g_rx[g_phase == 0][maske1]))
plt.grid()

plt.subplot(2,2,3)
maske2 = ~np.isnan(g_z2d[g_phase  == 1]) & ~np.isnan(g_rx[g_phase  == 1])
plt.hist2d(g_z2d[g_phase  == 1][maske2], g_rx[g_phase  == 1][maske2],bins=99)
plt.colorbar()
plt.title('mixed - corr: '+ corcor(g_z2d[g_phase  == 1][maske2], g_rx[g_phase  == 1][maske2]))
plt.grid()

plt.subplot(2,2,4)
maske3 = ~np.isnan(g_z2d[g_phase  == 2]) & ~np.isnan(g_rx[g_phase  == 2])
plt.hist2d(g_z2d[g_phase  == 2][maske3], g_rx[g_phase  == 2][maske3],bins=99)
plt.title('liquid - corr: '+corcor(g_z2d[g_phase  == 2][maske3], g_rx[g_phase  == 2][maske3]))
plt.xlabel('DPR')
plt.ylabel('radolan')
plt.grid()
plt.colorbar()
plt.show()


plt.figure(5)
plt.subplot(2,2,1)

plt.scatter(g_p2d[g_phase == 1], g_ry[g_phase == 1],
            alpha=0.2, color='blue', label='mixed')
plt.scatter(g_p2d[g_phase == 2], g_ry[g_phase == 2],
            alpha=0.2, color='red', label='liquid')
plt.scatter(g_p2d[g_phase == 0], g_ry[g_phase== 0],
            alpha=0.5, color='black', label='solid')
plt.xlabel('dpr')
plt.ylabel('radolan')

plt.legend()
plt.grid()

plt.subplot(2,2,2)
maske1 = ~np.isnan(g_p2d[g_phase == 0]) & ~np.isnan(g_ry[g_phase == 0])
plt.hist2d(g_p2d[g_phase  == 0][maske1], g_ry[g_phase == 0][maske1],bins=99)
plt.colorbar()
plt.title('solid - corr: '+ corcor(g_p2d[g_phase  == 0][maske1], g_ry[g_phase == 0][maske1]))
plt.grid()

plt.subplot(2,2,3)
maske2 = ~np.isnan(g_p2d[g_phase  == 1]) & ~np.isnan(g_ry[g_phase  == 1])
plt.hist2d(g_p2d[g_phase  == 1][maske2], g_ry[g_phase  == 1][maske2],bins=99)
plt.colorbar()
plt.title('mixed - corr: '+ corcor(g_p2d[g_phase  == 1][maske2], g_ry[g_phase  == 1][maske2]))
plt.grid()

plt.subplot(2,2,4)
maske3 = ~np.isnan(g_p2d[g_phase  == 2]) & ~np.isnan(g_ry[g_phase  == 2])
plt.hist2d(g_p2d[g_phase  == 2][maske3], g_ry[g_phase  == 2][maske3],bins=99)
plt.title('liquid - corr: '+corcor(g_p2d[g_phase  == 2][maske3], g_ry[g_phase  == 2][maske3]))
plt.xlabel('DPR')
plt.ylabel('radolan')
plt.grid()
plt.colorbar()
plt.show()


plt.figure(6)
A, B = g_z2d[g_bbh<1000], g_rx[g_bbh<1000]
plt.scatter(A, B, label=corcor(A,B) )
plt.legend()
plt.show()


plt.figure(7)
A, B = g_z2d[g_bbw>1000], g_rx[g_bbw>1000]
plt.scatter(A, B, label=corcor(A,B))
plt.legend()
plt.show()


for i in range(0,3100,100):
    print i
    A, B = g_z2d[(g_bbh>i) & (g_bbh <i+100)], g_rx[(g_bbh>i) & (g_bbh <i+100)]
    plt.scatter(A, B, label=str(i)+'- -'+corcor(A,B))
    plt.legend()

plt.show()