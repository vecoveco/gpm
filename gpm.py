'''Dieses Program soll dazu dienen die 
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Regenraten zu validieren.
Hier wird nur ein Ueberflug analysiert'''

#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------
'''Einlesen von Modulen'''
#--------------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob

'''ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/'''
#--------------------------------------------------------------------------------------------------------
### Pfad mit String ##
#--------------------------------------------------------------------------------------------------------
#Hier das Datum von den GPM Ueberfluegen nehmen
year = '2014'
m = '06'
d = '29'
#Hier die Uhrzeit der dazugehoeringen Boxpol RHI Scans
#ht:Stunde; mt:Minute; st:Sekunde
ht = '14'
mt = '50'
st = '00'



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


#--------------------------------------------------------------------------------------------------------
### Hier kann man gucken welchen Pfadtyp man dann zum einlesen benutzt, den oben oder den unten ##
#--------------------------------------------------------------------------------------------------------

#gpmrhi01_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_01/2014-06-29--14:54:52,00.mvol'
#gpmrhi02_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_02/2014-06-29--14:54:52,00.mvol'
#gpmrhi03_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_03/2014-06-29--14:54:52,00.mvol'
#gpmrhi04_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_04/2014-06-29--14:54:52,00.mvol'
#gpmrhi05_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_05/2014-06-29--14:54:52,00.mvol'
#gpmrhi06_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_06/2014-06-29--14:54:52,00.mvol'
#gpmrhi07_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_07/2014-06-29--14:54:52,00.mvol'
#gpmrhi08_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_08/2014-06-29--14:54:52,00.mvol'
#gpmrhi09_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_09/2014-06-29--14:54:52,00.mvol'
#gpmrhi10_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_10/2014-06-29--14:54:52,00.mvol'
#gpmrhi11_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_11/2014-06-29--14:54:52,00.mvol'
#gpmrhi12_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_12/2014-06-29--14:54:52,00.mvol'
#gpmrhi13_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_13/2014-06-29--14:54:52,00.mvol'
#gpmrhi14_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_14/2014-06-29--14:54:52,00.mvol'
#gpmrhi15_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_15/2014-06-29--14:54:52,00.mvol'

#gpmgmi_datapath = '/home/velibor/shkgpm/data/20140629/gprof/2A.GPM.GMI.GPROF2014v1-4.20140629-S140300-E153533.001897.V03C.HDF5'

#dpr_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.DPR.V5-20140827.20140629-S140300-E153533.001897.V03B.HDF5'
#ka_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.Ka.V5-20140829.20140629-S140300-E153533.001897.V03B.HDF5'
#ku_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.Ku.V5-20140829.20140629-S140300-E153533.001897.V03B.HDF5'
#dprgmi_datapath = '/home/velibor/shkgpm/data/20140629/radar/2B.GPM.DPRGMI.CORRA2014.20140629-S140300-E153533.001897.V03C.HDF5'


#ppi_datapath = "/automount/radar/scans/2014/2014-06/2014-06-29/ppi_1p5deg/2014-06-29--14:59:25,00.mvol"##
ppi_datapath=('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",00.mvol")

#--------------------------------------------------------------------------------------------------------
### Rhi BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

gpmrhi01 = h5py.File(pfad_boxpol_rhi01, 'r')


#--------------------------------------------------------------------------------------------------------
### PPI BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------
#ppi_datapath = "/automount/radar/scans/2014/2014-06/2014-06-29/ppi_1p5deg/2014-06-29--14:25:00,00.mvol"
#ppi_datapath = "/automount/radar/scans/2014/2014-06/2014-06-29/ppi_1p5deg/2014-06-29--14:59:25,00.mvol"##

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH = data['SCAN0']['ZH']['data']
r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']
#Umwandeln von Z in RR
Z = wradlib.trafo.idecibel(ZH)
R = wradlib.zr.z2r(Z, a=200., b=1.6) # Rainfall intensity using the Marshal-Palmer Z(R) parameters
depth = wradlib.trafo.r2depth(R, 300)# Convert to rainfall depth assuming a duration of 300 seconds
#--------------------------------------------------------------------------------------------------------
### gprof Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

gpmgmi = h5py.File(pfad_gprof_g, 'r')

#--------------------------------------------------------------------------------------------------------
### radar Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

dpr = h5py.File(pfad_radar_DPR, 'r')

ka = h5py.File(pfad_radar_Ka, 'r')
ku = h5py.File(pfad_radar_Ku, 'r')
dprgmi = h5py.File(pfad_radar_DPRGMI, 'r')

#--------------------------------------------------------------------------------------------------------
###---- GPROF ---- Einlesen von Oberflachenniederschlag und Lat/Lon von Gprof ##
#--------------------------------------------------------------------------------------------------------

gpmgmi.keys()
gpmgmi_S1=gpmgmi['S1']				
gprof_lat=gpmgmi_S1['Latitude']			#(2962, 221)
gprof_lon=gpmgmi_S1['Longitude']		#(2962, 221)
gprof_pp=gpmgmi_S1['surfacePrecipitation']   	#(2962, 221)

#--------------------------------------------------------------------------------------------------------
### ---- Rhi BoXPol ---- Einlesen... ##
#--------------------------------------------------------------------------------------------------------
'''
scan=gpmrhi01['scan0']
boxpol_zdr=scan['moment_9']   			#(58, 1800)

ppi1=ppi['scan0']
ppi_zdr=ppi1['moment_9']
ppi_uv=ppi1['moment_4']
'''
#--------------------------------------------------------------------------------------------------------
### ---- radar dpr ka ku dprgmi---- ##
#--------------------------------------------------------------------------------------------------------

dpr_NS = dpr['NS']				#Einlesen ueber HS und MS
dpr_SLV = dpr_NS['SLV']
dpr_pr = dpr_SLV['precipRate']			#(7934, 49, 176)
dpr_prs = dpr_SLV['precipRateNearSurface']	#(7934, 49)
dpr_lat = dpr_NS['Latitude']			#(7934, 49)
dpr_lon = dpr_NS['Longitude']			#(7934, 49)


#--------------------------------------------------------------------------------------------------------
### ------- In Arrays umwandeln ----------- ##
#--------------------------------------------------------------------------------------------------------

gprof_pp_a = np.array(gprof_pp)
gprof_lon_a = np.array(gprof_lon)
gprof_lat_a = np.array(gprof_lat)

dpr_lon_a = np.array(dpr_lon)
dpr_lat_a = np.array(dpr_lat)
dpr_pr_a = np.array(dpr_pr)
dpr_prs_a = np.array(dpr_prs)

#ppi_zdr_a=np.array(ppi_zdr)



#--------------------------------------------------------------------------------------------------------
### ------- Lon Lat Bestimmung ----------- ##
#--------------------------------------------------------------------------------------------------------
# Indexzugriff
#   gprof_lat_a[[Zeile][Zeilenelement]]
# Indexsuche
#   itemin= np.where((gprof_lat_a<58) & (gprof_lat_a>57))
bonn_lat1 = 49.9400
bonn_lat2 = 51.3500
bonn_lon1 = 6.40000
bonn_lon2 = 8.10000

ilat= np.where((gprof_lat_a>49.9400) & (gprof_lat_a<51.3500))
ilon= np.where((gprof_lon_a>6.40000) & (gprof_lon_a<8.10000))
lonstart = ilon[0][0]#erstes Element
lonend = ilon[0][-1]#letztes Element
latstart = ilat[0][0]
latend = ilat[0][-1]



gprof_pp_a[gprof_pp_a==-9999] = np.nan #NAN!

#--------------------------------------------------------------------------------------------------------
''' --------------------------------- Polar zu LonLat ------------------------------------- ''' 
#--------------------------------------------------------------------------------------------------------


radar_location = (lon_ppi, lat_ppi, alt_ppi) # (lon, lat, alt) in decimal degree and meters
elevation = 1.5 # in degree
azimuths = az # in degrees
ranges = r # in meters
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

gk3 = wradlib.georef.epsg_to_osr(31467)
x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
#ae = wradlib.georef.create_osr("aeqd", lon_0=radar_location[0], lat_0=radar_location[1])
#x, y = wradlib.georef.reproject(lon, lat, projection_target=ae)
xgrid, ygrid = wradlib.georef.reproject(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], projection_target=gk3)

#xgrid = np.linspace(x.min(), x.mean(), 100)
#ygrid = np.linspace(y.min(), y.mean(), 100)
#grid_xy = np.meshgrid(xgrid, ygrid)
#grid_xy = np.vstack((grid_xy[0].ravel(), grid_xy[1].ravel())).transpose()
grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), wradlib.ipol.Idw)#Linear, Idw, Nearest
gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)

#hier gridded umwandeln in RR

fig = plt.figure(figsize=(10,8))
ax = plt.subplot(111, aspect="equal")
pm = plt.pcolormesh(xgrid, ygrid, gridded)#xgrid, ygrid,
plt.xlim((x.min(),x.max()))
plt.ylim((y.min(),y.max()))
plt.colorbar(pm, shrink=0.75)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.savefig('interpoliert.png')
plt.close()

fig = plt.figure(figsize=(10,8))
ax = plt.subplot(111, aspect="equal")
pm = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], gridded)#xgrid, ygrid,
plt.xlim((lon.min(),lon.max()))
plt.ylim((lat.min(),lat.max()))
plt.colorbar(pm, shrink=0.75)
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
#plt.savefig('interpoliert2.png')
#plt.close()
plt.show()
#--------------------------------------------------------------------------------------------------------
''' --------------------------------- PLOTS ------------------------------------- ''' 
#--------------------------------------------------------------------------------------------------------


'''Vergleichsplot'''

#plt.figure(1)
fig = plt.figure(figsize=(13,10))
#Levels berechnen
maxv = np.max([np.max(gridded),np.max(np.ma.masked_invalid(gprof_pp_a)[latstart:latend])])
#
plt.subplot(221)

A = gridded
B = np.ma.masked_invalid(gprof_pp_a)[latstart:latend]
#Nullen entfernen
#A[A<1.0]=np.nan
#B[B<1.0]=np.nan
#Ausreiser entfernen
#A[A>np.nanstd(A)]=np.nan
#B[B>np.nanstd(B)]=np.nan

#Scatter mit regrssion
from scipy import stats
mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept
plt.plot(B,line,'r-',B,A,'o')
#plt.gca().set_aspect('equal', adjustable='box')		#control quadratic x y axis
#plt.axis('equal')						#control quadratic x y axis
plt.xlim(0,maxv)						#control quadratic x y axis
plt.ylim(0,maxv)						#control quadratic x y axis
plt.xlabel("gprof")
plt.ylabel("ppi BoxPol")
plt.title("Scatterplot (gprof/ppi), cor: " + str(r_value))


##
plt.subplot(222)
#rainrate
ax1, pm2 = wradlib.vis.plot_ppi(R,r,az,vmin=0,vmax=maxv)
cbar = plt.colorbar(pm2, shrink=0.75)
cbar.set_label("RR [mm/h]")
plt.xlim((-101000,101000))
plt.ylim((-101000,101000))
plt.xticks(())
plt.yticks(())
plt.xlabel("X Range [km]")
plt.ylabel("Y Range [km]")
plt.title(ppi_datapath[-28:-8])

plt.subplot(223)
#gprof rainrate
pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], np.ma.masked_invalid(gprof_pp_a)[latstart:latend],vmin=0,vmax=maxv)
plt.xlim((bonn_lon1,bonn_lon2))
plt.ylim((bonn_lat1,bonn_lat2))
plt.title(pfad_boxpol_rhi01[-28:-6])
cbar = plt.colorbar(pm2, shrink=.75)
cbar.set_label("GPROF surfacePrecipitation [mm/h]")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")

plt.subplot(224)
#ppi rainrate
pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], gridded,vmin=0,vmax=maxv)#xgrid, ygrid,
#plt.xlim((lon.min(),lon.max()))
#plt.ylim((lat.min(),lat.max()))
plt.xlim((bonn_lon1,bonn_lon2))
plt.ylim((bonn_lat1,bonn_lat2))
plt.title(ppi_datapath[-28:-8])
cbar = plt.colorbar(pm2, shrink=0.75)
cbar.set_label("Boxpol_ppi_interpolation RR [mm/h]")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")

#plt.savefig(ppi_datapath[-28:-8] + '_Vergleich.png')
#plt.close()
plt.show()




