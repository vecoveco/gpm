import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace

import matplotlib.pyplot as plt

para = ['RR_MS', 'RR_NS', 'RR_HS', 'REF_MS', 'REF_NS', 'REF_HS', 'SKILL_NS', 'SKILL_HS', 'SKILL_MS']

# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]

for ip in range(6):

    pp = ip


    pfad = ('/automount/ags/velibor/gpmdata/dumpdata/' + para[pp] + '/' + para[pp] + '*.hdf5')

    pfad_gpm = sorted(glob.glob(pfad))

    print 'Es sind ',len(pfad_gpm), 'Dateien vorhanden!'

    r_rad, r_sat = [], []
    rmax_rad, rmax_sat = [], []
    rmin_rad, rmin_sat = [], []
    rm_rad, rm_sat = [], []
    gx, gy = [], []
    rx, ry = [], []

    #len(pfad_gpm)

    for ii in range(len(pfad_gpm)):

        print ii,' von ', len(pfad_gpm)

        R = h5py.File(pfad_gpm[ii], 'r')
        rrr = np.array(R['dat_rad'])
        ggg = np.array(R['dat_sat'])
        x = np.array(R['x'])
        y = np.array(R['y'])
        rrr[rrr<=0] = np.nan
        ggg[ggg<=0] = np.nan

        rrr[rrr > 1000000000] = np.nan
        ggg[ggg > 1000000000] = np.nan

        rrr = rrr.reshape(rrr.shape[0]*rrr.shape[1])
        ggg = ggg.reshape(ggg.shape[0]*ggg.shape[1])
        x = x.reshape(x.shape[0] * x.shape[1])
        y = y.reshape(y.shape[0] * y.shape[1])

        r_rad = np.concatenate((r_rad, rrr), axis=0)
        r_sat = np.concatenate((r_sat, ggg), axis=0)
        gx = np.concatenate((gx, x), axis=0)
        gy = np.concatenate((gy, y), axis=0)


    print r_rad.shape, r_sat.shape, gx.shape, gy.shape
    np.save('/automount/ags/velibor/gpmdata/dumpdata/npy/'+str(para[pp]),(gx,gy,r_rad,r_sat))

g_bbh, g_bbw, g_pha, g_typ, g_x, g_y = [], [], [], [], [], []

for ii in range(len(pfad_gpm)):

    print ii,' von ', len(pfad_gpm)

    R = h5py.File(pfad_gpm[ii], 'r')
    bbh = np.array(R['sat_bbh'])
    bbw = np.array(R['sat_bbw'])
    pha = np.array(R['sat_phase'])
    #pp3 = np.array(R['sat_pp'])
    typ = np.array(R['sat_type'])
    x = np.array(R['x'])
    y = np.array(R['y'])

    # Daten reshapen
    bbh = bbh.reshape(bbh.shape[0]*bbh.shape[1])
    bbw = bbw.reshape(bbw.shape[0]*bbw.shape[1])
    pha = pha.reshape(pha.shape[0]*pha.shape[1])
    typ = typ.reshape(typ.shape[0]*typ.shape[1])
    x = x.reshape(x.shape[0] * x.shape[1])
    y = y.reshape(y.shape[0] * y.shape[1])

    #Daten Conecten
    g_bbh = np.concatenate((g_bbh, bbh), axis=0)
    g_bbw = np.concatenate((g_bbw, bbw), axis=0)
    g_pha = np.concatenate((g_pha, pha), axis=0)
    g_typ = np.concatenate((g_typ, typ), axis=0)
    g_x = np.concatenate((g_x, x), axis=0)
    g_y = np.concatenate((g_y, y), axis=0)


print g_bbw.shape, g_bbh.shape, g_x.shape, g_y.shape, g_pha.shape, g_typ.shape