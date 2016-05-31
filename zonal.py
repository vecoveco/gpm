# Programm um boxpol auf gprof Gitter zu proizieren
# mit zonal statistics.
# Validierung von Regen Raten beider Produkte

import os
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import from_levels_and_colors
import matplotlib.patches as patches
import datetime as dt
import matplotlib.pyplot  as plt
import wradlib
import glob
import h5py
from osgeo import osr

print "...zonal.py wurde gestartet..."

LISTE = ("20150128172208", "20140629145000", "20140629145925", "20140921070500", "20140921071058",
         "20141007023744")
LISTE = sorted(LISTE)

ZP = LISTE[0]
year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]

pfad = ('/home/velibor/shkgpm/data/' + year + m + d + '/radar/*.HDF5')
pfad_radar = sorted(glob.glob(pfad))
pfad_radar_Ku = pfad_radar[2]
pfad_radar_DPR = pfad_radar[0]
pfad_radar_DPRGMI = pfad_radar[3]
pfad_radar_Ka = pfad_radar[1]

pfad1 = ('/home/velibor/shkgpm/data/' + year + m + d + '/boxpol/gpm_rhi_01/*.mvol')
pfad_boxpol = glob.glob(pfad1)
pfad_boxpol_rhi01 = pfad_boxpol[0]

pfad2 = ('/home/velibor/shkgpm/data/' + year + m + d + '/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

ppi_datapath = (
'/automount/radar-archiv/scans/' + year + "/" + year + "-" + m + "/" + year + "-" + m + "-" + d + "/ppi_1p5deg/" + year + "-" + m + "-" + d + "--" + ht + ":" + mt + ":" + st + ",00.mvol")

# Rhi BoxPol Daten einlesen
gpmrhi01 = h5py.File(pfad_boxpol_rhi01, 'r')

# PPI BoxPol Daten einlesen
ppi = h5py.File(ppi_datapath, 'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

# setup OSR objects
proj_gk = osr.SpatialReference()
proj_gk.ImportFromEPSG(31466)
proj_ll = osr.SpatialReference()
proj_ll.ImportFromEPSG(4326)

ZH = data['SCAN0']['ZH']['data']
PHIDP = data['SCAN0']['PHIDP']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']
# Umwandeln von Z in RR Marshal-Palmer Z(R)
Z = wradlib.trafo.idecibel(ZH)
R = wradlib.zr.z2r(Z, a=200., b=1.6)


rays = az.shape[0]
bins = r.shape[0]

# create polar grid polygon vertices in lat,lon
radar_ll = wradlib.georef.polar2polyvert(r, az, (lon_ppi, lat_ppi))
# project ll grids to GK2
radar_gk = wradlib.georef.reproject(radar_ll, projection_source=proj_ll,
                                    projection_target=proj_gk)
# reshape
radar_gk.shape = (rays, bins, 5, 2)

### gprof Daten einlesen
gpmgmi = h5py.File(pfad_gprof_g, 'r')

# GPROF ---- Einlesen von Oberflachenniederschlag und Lat/Lon von Gprof
gpmgmi.keys()
gpmgmi_S1 = gpmgmi['S1']
gprof_lat = gpmgmi_S1['Latitude']
gprof_lon = gpmgmi_S1['Longitude']
gprof_pp = gpmgmi_S1['surfacePrecipitation']

gprof_pp_a = np.array(gprof_pp)
gprof_lon_a = np.array(gprof_lon)
gprof_lat_a = np.array(gprof_lat)
print("gprof_lon_a.shape", gprof_lat_a.shape)

ilat= np.where((gprof_lat_a>49.9400) & (gprof_lat_a<51.3500))
ilon= np.where((gprof_lon_a>6.40000) & (gprof_lon_a<8.10000))
lonstart = ilon[0][0]	#erstes Element
lonend = ilon[0][-1]	#letztes Element
latstart = ilat[0][0]
latend = ilat[0][-1]
gp_lon = gprof_lon_a[latstart:latend]
ilon= np.where((gp_lon>6.40000) & (gp_lon<8.10000)) #Hier eventuell doppelt ? s.o.
lonstart = ilon[0][0]	#erstes Element
lonend = ilon[0][-1]	#letztes Element
gp_lon1 = gp_lon[lonstart:lonend] # Index Eingrenzung HIER KOORDINATEN RAUSSUCHEN!

xgrid, ygrid = wradlib.georef.reproject(gprof_lon_a[latstart:latend],
                                        gprof_lat_a[latstart:latend], projection_target=proj_gk)

gprof_pp_a[gprof_pp_a == -9999] = np.nan

gprof_git = gprof_pp_a[latstart:latend]
gprof_git1 = gprof_lon[lonstart:lonend]

gprof_lon = gprof_lon_a[latstart:latend]
gprof_lon1 = gprof_lon[lonstart:lonend]  # Werte von gprof eingegrenzt mit lon lat Limits

gprof_lat = gprof_lat_a[latstart:latend]
gprof_lat1 = gprof_lat[lonstart:lonend]
# Nur Ideen
print ('SHAPE: gprof_lon1')
print (gprof_lon1.shape, gprof_lon1.dtype)


#Todo: berechnetes Grid Fehlerhaft
gprof_gitter_lon = np.empty((98,221,5))
gprof_gitter_lat = np.empty((98,221,5))
for i in range(0,98,1):
    for j in range(0,221,1):
        if j < 220 and i < 97:
            gprof_gitter_lon [[i],[j],[0]] = ((gprof_lon1[[i+1],[j]] + (gprof_lon1[[i+1],[j + 1]] - gprof_lon1[[i],[j]])/2) - (gprof_lon1[[i],[j]] + (gprof_lon1[[i],[j + 1]] - gprof_lon1[[i],[j]])/2))/2
            gprof_gitter_lat [[i],[j],[0]] = ((gprof_lat1[[i+1],[j]] + (gprof_lat1[[i+1],[j + 1]] - gprof_lat1[[i],[j]])/2) - (gprof_lat1[[i],[j]] + (gprof_lat1[[i],[j + 1]] - gprof_lat1[[i],[j]])/2))/2
            # erste Koo
            gprof_gitter_lon [[i],[j],[1]] = ((gprof_lon1[[i],[j-1]] + (gprof_lon1[[i - 1],[j-1]] - gprof_lon1[[i],[j-1]])/2) - (gprof_lon1[[i],[j]] + (gprof_lon1[[i + 1],[j]] - gprof_lon1[[i],[j]])/2))/2
            gprof_gitter_lat [[i],[j],[1]] = ((gprof_lat1[[i],[j-1]] + (gprof_lat1[[i - 1],[j-1]] - gprof_lat1[[i],[j-1]])/2) - (gprof_lon1[[i],[j]] + (gprof_lon1[[i + 1],[j]] - gprof_lon1[[i],[j]])/2))/2
            # zweite Koo
            gprof_gitter_lon [[i],[j],[2]] = ((gprof_lon1[[i-1],[j]] - (gprof_lon1[[i-1],[j - 1]] - gprof_lon1[[i-1],[j]])/2) - (gprof_lon1[[i],[j]] - (gprof_lon1[[i],[j + 1]] - gprof_lon1[[i],[j]])/2))/2
            gprof_gitter_lat [[i],[j],[2]] = ((gprof_lat1[[i-1],[j]] - (gprof_lat1[[i-1],[j - 1]] - gprof_lat1[[i-1],[j]])/2) - (gprof_lat1[[i],[j]] - (gprof_lat1[[i],[j + 1]] - gprof_lat1[[i],[j]])/2))/2
            # dritte Koo
            gprof_gitter_lon [[i],[j],[3]] = ((gprof_lon1[[i],[j+1]] - (gprof_lon1[[i + 1],[j+1]] - gprof_lon1[[i],[j+1]])/2) - (gprof_lon1[[i],[j]] - (gprof_lon1[[i + 1],[j]] - gprof_lon1[[i],[j]])/2))/2
            gprof_gitter_lat [[i],[j],[3]] = ((gprof_lat1[[i],[j+1]] - (gprof_lat1[[i + 1],[j+1]] - gprof_lat1[[i],[j+1]])/2) - (gprof_lat1[[i],[j]] - (gprof_lat1[[i + 1],[j]] - gprof_lat1[[i],[j]])/2))/2
            # vierte Koo
            gprof_gitter_lon [[i],[j],[0]] = ((gprof_lon1[[i+1],[j]] + (gprof_lon1[[i+1],[j + 1]] - gprof_lon1[[i],[j]])/2) - (gprof_lon1[[i],[j]] + (gprof_lon1[[i],[j + 1]] - gprof_lon1[[i],[j]])/2))/2
            gprof_gitter_lat [[i],[j],[0]] = ((gprof_lat1[[i+1],[j]] + (gprof_lat1[[i+1],[j + 1]] - gprof_lat1[[i],[j]])/2) - (gprof_lat1[[i],[j]] + (gprof_lat1[[i],[j + 1]] - gprof_lat1[[i],[j]])/2))/2
            # erste Koo
        else:   # Randbedinungungen: bei j=221 und i=98 gibt es kein j+1 und i +1 Idee ji davor
            #Todo: Randbedingungen bei j = 0 und i= 0 auch neu erstellen !
            gprof_gitter_lon [[i],[j],[0]] = gprof_lon1[[i],[j]] + (gprof_lon1[[i],[j - 1]] - gprof_lon1[[i],[j]])/2
            gprof_gitter_lat [[i],[j],[0]] = gprof_lat1[[i],[j]] + (gprof_lat1[[i],[j - 1]] - gprof_lat1[[i],[j]])/2
            # erste Koo bei RB
            gprof_gitter_lon [[i],[j],[1]] = gprof_lon1[[i],[j]] + (gprof_lon1[[i - 1],[j]] - gprof_lon1[[i],[j]])/2
            gprof_gitter_lat [[i],[j],[1]] = gprof_lat1[[i],[j]] + (gprof_lat1[[i - 1],[j]] - gprof_lat1[[i],[j]])/2
            # zweite Koo bei RB
            gprof_gitter_lon [[i],[j],[2]] = gprof_lon1[[i],[j]] - (gprof_lon1[[i],[j - 1]] - gprof_lon1[[i],[j]])/2
            gprof_gitter_lat [[i],[j],[2]] = gprof_lat1[[i],[j]] - (gprof_lat1[[i],[j - 1]] - gprof_lat1[[i],[j]])/2
            # dritte Koo bei RB
            gprof_gitter_lon [[i],[j],[3]] = gprof_lon1[[i],[j]] - (gprof_lon1[[i - 1],[j]] - gprof_lon1[[i],[j]])/2
            gprof_gitter_lat [[i],[j],[3]] = gprof_lat1[[i],[j]] - (gprof_lat1[[i - 1],[j]] - gprof_lat1[[i],[j]])/2
            # vierte Koo bei RB
            gprof_gitter_lon [[i],[j],[4]] = gprof_lon1[[i],[j]] + (gprof_lon1[[i],[j - 1]] - gprof_lon1[[i],[j]])/2
            gprof_gitter_lat [[i],[j],[4]] = gprof_lat1[[i],[j]] + (gprof_lat1[[i],[j - 1]] - gprof_lat1[[i],[j]])/2
            # erste Koo bei RB

gprof_gitter_lon = gprof_gitter_lon[...,np.newaxis]
gprof_gitter_lat = gprof_gitter_lat[...,np.newaxis]
print(gprof_gitter_lon.shape)
gprof_gitter = np.concatenate((gprof_gitter_lon, gprof_gitter_lat ), axis=3)
print(gprof_gitter.shape)
gprof_gitter = gprof_gitter.reshape(98*221,5,2)
print(gprof_gitter.shape)
print(gprof_gitter[0,0,:])


zd = wradlib.zonalstats.ZonalDataPoly(radar_ll, gprof_gitter[::100], srs=proj_ll, buf=0.)
zd.dump_vector('gprof_isect')

gprof_raw_gitter = np.dstack((gprof_lon_a, gprof_lat_a))
print(gprof_raw_gitter.shape)
zd1 = wradlib.zonalstats.ZonalDataPoint(gprof_raw_gitter.reshape(2962*221,2), gprof_gitter[:10], srs=proj_ll, buf=0.)
zd1.dump_vector('gprof_isect1')

print ('gitter100  ',gprof_gitter.shape)
#Todo: BSP PLOT!
#plt.pcolormesh(gprof_gitter[:100])
#plt.show()