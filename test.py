
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

print ("...test.py wurde gestartet...")

LISTE = ("20140629145000", "20140629145925", "20140921070500", "20140921071058",
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

print ("...Radar und Satellitdaten wurden eingelesen...")


# ============================= EINLESEN DER SHP DATEN ======================================================= #

zd = wradlib.zonalstats.ZonalDataPoly('/user/velibor/SHKGPM/data/test_zonal_poly_cart')
# isecs = zd.get_isecs
#z = wradlib.zonalstats.ZonalDataPoly.load_vector('/user/velibor/SHKGPM/data/test_zonal_poly_cart/trg.shp')
print("zd loaded")
obj3 = wradlib.zonalstats.GridCellsToPoly(zd) #Dauert sehr lange
obj3.zdata.dst.set_attribute('ix', obj3.ix)
obj3.zdata.dst.set_attribute('w', obj3.w)

print(dir(obj3))
obj3.zdata.dump_vector('/user/velibor/SHKGPM/data/test_zonal_poly')
#print(zd.dst.data.shape)#.shape, zd.trg.shape, zd.dst.shape)

avg1 = obj3.mean(R.ravel())
print(avg1)

#exit(9)

'''Dauer!?
# get intersections as numpy array
isecs = zd.isecs

print(isecs.shape)
'''

''' Dauer!?
gc = wradlib.zonalstats.GridCellsToPoly(zd)
count = radar_gk.shape[0]
data = 1000000. / np.array(range(count))
# calculate mean and variance
mean = gc.mean(data)
var = gc.var(data)

print("Average:", mean)
print("Variance:", var)
'''
"""
# GProf Gitter ueber Radar
dst_shp = wradlib.util.get_wradlib_data_file('/user/velibor/SHKGPM/data/test_zonal_poly_cart/dst.shp')
d_dataset, d_inLayer = wradlib.io.open_shape(dst_shp)
d_cats, d_keys = wradlib.georef.get_shape_coordinates(d_inLayer)

dbox = d_inLayer.GetExtent()


# Radar Gitter
src_shp = wradlib.util.get_wradlib_data_file('/user/velibor/SHKGPM/data/test_zonal_poly_cart/src.shp')
s_dataset, s_inLayer = wradlib.io.open_shape(src_shp)
s_cats, s_keys = wradlib.georef.get_shape_coordinates(s_inLayer)

sbox = s_inLayer.GetExtent()

# Gprof ganzes  GitterPoly
trg_shp = wradlib.util.get_wradlib_data_file('/user/velibor/SHKGPM/data/test_zonal_poly_cart/trg.shp')
t_dataset, t_inLayer = wradlib.io.open_shape(trg_shp)
t_cats, t_keys = wradlib.georef.get_shape_coordinates(t_inLayer)

tbox = t_inLayer.GetExtent()

print(s_dataset)
print('d_cats: ',d_cats.shape,' d_keys: ', d_keys, ' d_inLayer: ', d_inLayer)
print('s_cats: ',s_cats.shape,' s_keys: ', s_keys)
print('t_cats: ',t_cats.shape,' t_keys: ', t_keys)
print('sbox: ',sbox)
print('dbox: ',dbox)
print('tbox: ',dbox)
"""





############# PLOT #######################
fig = plt.figure(figsize=(13,10))
plt.subplot(211)
pm2 = wradlib.vis.plot_ppi(R,r,az)
bbox = wradlib.zonalstats.get_bbox(t_cats[1][:, 1], t_cats[1][:, 1])
#plt.xlim((bonn_lon1,bonn_lon2))
#plt.ylim((bonn_lat1,bonn_lat2))
plt.subplot(212)
#pm2 = plt.pcolormesh(np.ma.masked_invalid(gprof_pp_a[sbox]))
#plt.savefig('/user/velibor/SHKGPM/data/plot/' + ppi_datapath[-28:-8] + '_Vergleich.png')
#plt.close()
i=1
bbox = wradlib.zonalstats.get_bbox(t_cats, t_cats)
plt.show()

