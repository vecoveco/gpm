import h5py
import numpy as np
import wradlib
import glob
from scipy import stats, linspace

import matplotlib.pyplot as plt

para = ['RR_MS', 'RR_NS', 'RR_HS', 'REF_MS', 'REF_NS', 'REF_HS']

# Threshold ATBD GPM 2016
thresh = [0.2, 0.5, 0.2, 12.0, 18.0, 12.0]

pp = 4

#TH = thresh[pp]
#TH = 0.15
#pfad = ('/automount/ags/velibor/gpmdata/dumpdata/RR_MS/RR_MS*.hdf5')
pfad = ('/automount/ags/velibor/gpmdata/dumpdata/' + para[pp] + '/' + para[pp] + '*.hdf5')

pfad_gpm = sorted(glob.glob(pfad))

print 'Es sind ',len(pfad_gpm), 'Dateien vorhanden!'

sommer, winter, = ['07','08','06'], ['12','01','02']
r_rad, r_sat = [], []
rmax_rad, rmax_sat = [], []
rmin_rad, rmin_sat = [], []
rm_rad, rm_sat = [], []
gx, gy = [], []
rx, ry = [], []

print len(pfad_gpm)

for ii in range(len(pfad_gpm)):
    print ii


    #if pfad_gpm[ii][56:58] in sommer:
    R = h5py.File(pfad_gpm[ii], 'r')
    rrr = np.array(R['dat_rad'])
    ggg = np.array(R['dat_sat'])
    x=  np.array(R['x'])
    y = np.array(R['y'])
    rrr[rrr<=0] = np.nan
    ggg[ggg<=0] = np.nan

    rrr[rrr>1000000000] = np.nan
    ggg[ggg>1000000000] = np.nan

    rrr = rrr.reshape(rrr.shape[0]*rrr.shape[1])
    ggg = ggg.reshape(ggg.shape[0]*ggg.shape[1])
    x = x.reshape(x.shape[0]*x.shape[1])
    y = y.reshape(y.shape[0]*y.shape[1])
    #print rrr.shape


    r_rad = np.concatenate((r_rad, rrr), axis=0)
    r_sat = np.concatenate((r_sat, ggg), axis=0)

    gx = np.concatenate((gx, x), axis=0)
    gy = np.concatenate((gy, y), axis=0)

    rmax_rad.append(np.nanmax(rrr))
    rmax_sat.append(np.nanmax(ggg))

    rmin_rad.append(np.nanmin(rrr))
    rmin_sat.append(np.nanmin(ggg))

    rm_rad.append(np.nanmean(rrr))
    rm_sat.append(np.nanmean(ggg))
    #else:
    #    pass
    #print r_rad.shape

thh = np.nanmax([np.nanmin(r_rad), np.nanmin(r_sat)])
print 'Min of sat and rad -------------->>>>>  ',thh
TH = 10

r_rad[r_rad<=TH] = np.nan
r_sat[r_sat<=TH] = np.nan

rad_coast = r_rad[np.where(gy >= -4000)]
sat_coast = r_sat[np.where(gy >= -4000)]

rad_complex = r_rad[np.where(gy <= -4000)]
sat_complex = r_sat[np.where(gy <= -4000)]

from satlib import validation_plot_log
from satlib import validation_plot
validation_plot(r_sat, r_rad, TH)
plt.title(str(para[pp])+ '- TH: ' + str(TH))
plt.show()
#plt.savefig('/automount/ags/velibor/plot/validation/validation'+para[pp]+'.png' )
#plt.close()
validation_plot(sat_coast, rad_coast, TH)
plt.title('Coast:_'+str(para[pp])+ '- TH: ' + str(TH))
plt.show()

validation_plot(sat_complex, rad_complex, TH)
plt.title('Complex:_'+str(para[pp])+ '- TH: ' + str(TH))
plt.show()


from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders
from pcc import plot_radar


fig = plt.figure(figsize=(12,12))
ax1 = fig.add_subplot(111, aspect='equal')
plt.hist2d(gx[np.where(r_rad > TH)],gy[np.where(r_rad > TH)],bins=100)
cbar = plt.colorbar()
cbar.set_label('# of Overpases with TH over'+ str(TH), rotation=270)
plt.grid()
plt.title('GPM NS REF Hist')
plot_borders(ax1)
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']

plot_radar(bonnlon, bonnlat, ax1, reproject=True, cband=False,col='black')

plt.show()
