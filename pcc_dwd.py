"""

Einlesen von dwd C-band Radardaten

"""


import h5py
import numpy as np
import wradlib
import glob
import wradlib as wrl
from osgeo import osr
from pcc import get_time_of_gpm
from pcc import cut_the_swath
import matplotlib.pyplot as plt

pfad_dwd_nhb = glob.glob('/automount/radar/dwd/dx/2014/2014-10/2014-10-07/raa00-dx_*-1410070235-nhb---bin')[0]
pfad_dwd_oft = glob.glob('/automount/radar/dwd/dx/2014/2014-10/2014-10-07/raa00-dx_*-1410070235-oft---bin')[0]
pfad_dwd_ess = glob.glob('/automount/radar/dwd/dx/2014/2014-10/2014-10-07/raa00-dx_*-1410070235-ess---bin')[0]



one_scan_nhb, attributes = wradlib.io.readDX(pfad_dwd_nhb)
one_scan_oft, attributes = wradlib.io.readDX(pfad_dwd_oft)
one_scan_ess, attributes = wradlib.io.readDX(pfad_dwd_ess)


pfad_boxpol = '/automount/radar-archiv/scans/2014/2014-10/2014-10-07/ppi_1p5deg/2014-10-07--02:37:44,00.mvol'
pfad_juxpol = '/automount/radar-archiv/scans_juelich/2014/2014-10/2014-10-07/DWD_Vol_2/2014-10-07--02:35:00,00.mvol'

data, attrs = wradlib.io.read_GAMIC_hdf5(pfad_boxpol)

ZH = data['SCAN0']['ZH']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']

radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)



dataj, attrsj = wradlib.io.read_GAMIC_hdf5(pfad_juxpol)

ZHj = dataj['SCAN0']['ZH']['data']
rj = attrsj['SCAN0']['r']
azj = attrsj['SCAN0']['az']
lon_ppij = attrsj['VOL']['Longitude']
lat_ppij = attrsj['VOL']['Latitude']
alt_ppij = attrsj['VOL']['Height']

radar_locationj = (lon_ppij, lat_ppij, alt_ppij)
elevationj = 1.5
azimuthsj = azj
rangesj = rj
polargridj = np.meshgrid(rangesj, azimuthsj)
lonj, latj, altj = wradlib.georef.polar2lonlatalt_n(polargridj[0], polargridj[1], elevationj, radar_locationj)


plt.subplot(2,2,2)
ax, pm = wradlib.vis.plot_ppi(one_scan_ess, vmin=0, vmax=50)
# add a colorbar with label
cbar = plt.colorbar(pm, shrink=0.75)
cbar.set_label("Reflectivity (dBZ)")
plt.title('Essen')

plt.subplot(2,2,3)
ax, pm = wradlib.vis.plot_ppi(one_scan_nhb, vmin=0, vmax=50)
# add a colorbar with label
cbar = plt.colorbar(pm, shrink=0.75)
cbar.set_label("Reflectivity (dBZ)")
plt.title('Neuheilenbach')


plt.subplot(2,2,4)
ax, pm = wradlib.vis.plot_ppi(one_scan_oft, vmin=0, vmax=50)
# add a colorbar with label
cbar = plt.colorbar(pm, shrink=0.75)
cbar.set_label("Reflectivity (dBZ)")
plt.title('Offenthal')

plt.show()


from pcc import get_radar_locations

radars = get_radar_locations()

ess_lon, ess_lat, ess_alt = radars['ESS']['lon'], radars['ESS']['lat'], radars['ESS']['alt']
oft_lon, oft_lat, oft_alt = radars['OFT']['lon'], radars['OFT']['lat'], radars['OFT']['alt']
nhb_lon, nhb_lat, nhb_alt = radars['NHB']['lon'], radars['NHB']['lat'], radars['NHB']['alt']


radar_location_oft = (oft_lon, oft_lat, oft_alt)
elevation_oft = 1.5
azimuths_oft = attributes['azim']
ranges_oft = np.arange(0, one_scan_oft.shape[1]) * 1000
polargrid = np.meshgrid(ranges_oft, azimuths_oft)
lonoft, latoft, altoft = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation_oft, radar_location_oft)

radar_location_nhb = (nhb_lon, nhb_lat, nhb_alt)
elevation_nhb = 1.5
azimuths_nhb = attributes['azim']
ranges_nhb = np.arange(0, one_scan_oft.shape[1]) * 1000
polargrid = np.meshgrid(ranges_nhb, azimuths_nhb)
lonnhb, latnhb, altnhb = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation_nhb, radar_location_nhb)

radar_location_ess = (ess_lon, ess_lat, ess_alt)
elevation_ess = 1.5
azimuths_ess = attributes['azim']
ranges_ess = np.arange(0, one_scan_oft.shape[1]) * 1000
polargrid = np.meshgrid(ranges_ess, azimuths_ess)
loness, latess, altess = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation_ess, radar_location_ess)


from pcc import get_miub_cmap
plt.subplot(2,2,1)
plt.pcolormesh(loness, latess, one_scan_ess, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(ess_lon, ess_lat)

plt.pcolormesh(lonnhb, latnhb, one_scan_nhb, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(nhb_lon, nhb_lat)

plt.pcolormesh(lonoft, latoft, one_scan_oft, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(oft_lon, oft_lat)

plt.subplot(2,2,2)

plt.pcolormesh(lonnhb, latnhb, one_scan_nhb, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(nhb_lon, nhb_lat)

plt.pcolormesh(lonoft, latoft, one_scan_oft, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(oft_lon, oft_lat)

plt.pcolormesh(loness, latess, one_scan_ess, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(ess_lon, ess_lat)

plt.subplot(2,2,3)
plt.pcolormesh(lonoft, latoft, one_scan_oft, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(oft_lon, oft_lat)

plt.pcolormesh(loness, latess, one_scan_ess, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(ess_lon, ess_lat)

plt.pcolormesh(lonnhb, latnhb, one_scan_nhb, vmin=0, vmax=50)#, cmap=get_miub_cmap())
plt.scatter(nhb_lon, nhb_lat)


#plt.pcolormesh(lon, lat, ZH, vmin=0, vmax=50, cmap=get_miub_cmap())
#plt.scatter(lon_ppi, lat_ppi)

#plt.pcolormesh(lonj, latj, ZHj, vmin=0, vmax=50, cmap=get_miub_cmap())
#plt.scatter(lon_ppij, lat_ppij)

plt.grid()
plt.show()