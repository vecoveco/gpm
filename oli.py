'''auslesen der hdf5 daten'''

#!/usr/bin/python

#--------------------------------------------------------------------------------------------------------
### module ##
#--------------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import csv
import matplotlib.pyplot as plt
import glob

np.set_printoptions(precision=4)

ZP = '20140921'; gpm_time = '2014-09-21'

year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
ye = ZP[2:4]

## GPM DPR
## ----------------------------

datei_gprof = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
#--------------------------------------------------------------------------------------------------------
### pfad ##
#--------------------------------------------------------------------------------------------------------
pfad_gprof = glob.glob(datei_gprof)
pfad_gprof_g = pfad_gprof[0]

#--------------------------------------------------------------------------------------------------------
### gprof-daten einlesen ##
#--------------------------------------------------------------------------------------------------------

gprof_dat = h5py.File(pfad_gprof_g, 'r')

#--------------------------------------------------------------------------------------------------------
### auswaehlen von lon lat und sp aus ordner S1 ##
#--------------------------------------------------------------------------------------------------------

gprof_dat.keys()
gpm_dat_S1=gprof_dat['S1']
gprof_lat=gpm_dat_S1['Latitude']		#(2962, 221)
gprof_lon=gpm_dat_S1['Longitude']		#(2962, 221)
gprof_sp=gpm_dat_S1['surfacePrecipitation']   	#(2962, 221)

#--------------------------------------------------------------------------------------------------------
### ------- in arrays umwandeln ----------- ##
#--------------------------------------------------------------------------------------------------------
gprof_lon_a = np.array(gprof_lon)
gprof_lat_a = np.array(gprof_lat)

gprof_sp_a = np.array(gprof_sp)

print gprof_sp_a.shape
#--------------------------------------------------------------------------------------------------------
### ------- lon lat definieren und auswerten, sowie "messfehler" aussortieren ----------- ##
#--------------------------------------------------------------------------------------------------------


# definition der flaeche
latbn1 = 49.5392
latbn2 = 50.8689
lonbn1 = 6.94964
lonbn2 = 7.30207

#wippen Koordinaten
wippe_lat = np.array([50.7455,50.7294,50.7114,50.6886,50.7689,50.6911,50.6696,50.6630,
             50.6815,50.6570,50.7118,50.7369,50.7502,50.7517,50.7269,50.7238,
             50.6392,50.7049,50.7269,50.7092,50.7127,50.6773])
wippe_lon = np.array([7.07457,7.08032,7.12037,7.08439,7.06215,7.15079,7.18338,7.14733,
             7.13585,7.12403, 7.14917,7.12871,7.20207,7.16625,7.18750,
             7.14808,7.12054,7.17520,7.04964,7.07456,7.06453,7.06609])

# auswerten
ilat = np.where((gprof_lat_a > latbn1) & (gprof_lat_a < latbn2))
ilon = np.where((gprof_lon_a > lonbn1) & (gprof_lon_a < lonbn2))

# loeschen der messfehler
gprof_sp_a[gprof_sp_a==-9999] = 0    #np.nan


# elementzuweisung
latstart = ilat[0][0]
latend = ilat[0][-1]
lonstart = ilon[0][0]
lonend = ilon[0][-1]



#--------------------------------------------------------------------------------------------------------
### ------- plot ----------- ##
#--------------------------------------------------------------------------------------------------------

pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], np.ma.masked_invalid(gprof_sp_a)[latstart:latend],vmin=0,vmax=0.3)
plt.xlim((lonbn1,lonbn2))
plt.ylim((latbn1,latbn2))
plt.title('GPROF GPM')
cbar = plt.colorbar(pm2, shrink=.75)
cbar.set_label('SurfacePrecipitation in [mm/h]')
plt.scatter(wippe_lon, wippe_lat)
plt.xlabel('Ost')
plt.ylabel('Nord')
plt.show()
