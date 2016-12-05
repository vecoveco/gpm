"""

Program zu erkennen der Position und Zeit des GPM Satelliten im Orbit

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib as mpl
import wradlib as wrl
from osgeo import osr
import os
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Zeitstempel nach YYYYMMDDhhmmss
# 20140921070500, 20150128170000 "20140629145000","20140629145925","20141007023744"
ZP = '20141007023744'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read GPROF
## ------------
#pfad2 = ('/home/velibor/shkgpm/data/20140921/gprof/*.HDF5')
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')

pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

print pfad_gprof_g

gpmgmi = h5py.File(pfad_gprof_g, 'r')

gpmgmi.keys()
gpmgmi_S1=gpmgmi['S1']
gprof_lat=np.array(gpmgmi_S1['Latitude'])
gprof_lon=np.array(gpmgmi_S1['Longitude'])
gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
gprof_pp[gprof_pp<=0] = np.nan


SC_lat =  np.array(gpmgmi_S1['SCstatus']['SClatitude'])
SC_lon =  np.array(gpmgmi_S1['SCstatus']['SClongitude'])
SC_alt =  np.array(gpmgmi_S1['SCstatus']['SCaltitude'])  # km
SC_ori =  np.array(gpmgmi_S1['SCstatus']['SCorientation'])

sc_dom =  np.array(gpmgmi_S1['ScanTime']['DayOfMonth'])
sc_doy =  np.array(gpmgmi_S1['ScanTime']['DayOfYear'])
sc_h =  np.array(gpmgmi_S1['ScanTime']['Hour'])
sc_ms =  np.array(gpmgmi_S1['ScanTime']['MilliSecond'])
sc_m =  np.array(gpmgmi_S1['ScanTime']['Minute'])
sc_M =  np.array(gpmgmi_S1['ScanTime']['Month'])
sc_sec =  np.array(gpmgmi_S1['ScanTime']['Second'])
sc_sod =  np.array(gpmgmi_S1['ScanTime']['SecondOfDay'])
sc_y =  np.array(gpmgmi_S1['ScanTime']['Year'])


#### Position fur Deutschland erstellen
bonn_lat1 = 47.00
bonn_lat2 = 55.00
bonn_lon1 = 6.00
bonn_lon2 = 15.00

ilat= np.where((gprof_lat>bonn_lat1) & (gprof_lat<bonn_lat2))
ilon= np.where((gprof_lon>bonn_lon1) & (gprof_lon<bonn_lon2))

lonstart = ilon[0][0]	#erstes Element
lonend = ilon[0][-1]	#letztes Element
latstart = ilat[0][0]
latend = ilat[0][-1]

#if SC_ori[0]==180:
#    offset = -25
#else:
#    offset = 50

pos = ((lonstart + lonend) /2) #+ offset

#Mittelpunkt von GMI Path
lon_scan = gprof_lon[pos,gprof_lon.shape[1]/2]
lat_scan = gprof_lat[pos,gprof_lat.shape[1]/2]

grad = 0.1
np.where((lon_scan+grad>=SC_lon)&(lon_scan-grad<=SC_lon))

dx, dy = [47,55,55,47], [6,6,15,15]


#### PLOT ####
for i in [1,2,3,4]:
    fig1 = plt.figure()
    pos = pos + i
    plt.plot(SC_lon, SC_lat, label='Satellitenbahn')
    plt.scatter(SC_lon[pos], SC_lat[pos], color='magenta', label='Satellitenposition')

    #plt.scatter((gprof_lon[pos,-1]+gprof_lon[pos,0])/2, (gprof_lat[pos,-1]+gprof_lat[pos,0])/2, label= 'Position GMI-Scan')

    plt.scatter(gprof_lon[pos,gprof_lon.shape[1]/2],
                gprof_lat[pos,gprof_lat.shape[1]/2],
                color='black', lw = 4, label='GMI Scan Mittel',
                zorder=2)


    plt.scatter(gprof_lon[pos,:], gprof_lat[pos,:], color='red', lw = 0.1, label='GMI Scan')
    #plt.scatter(gprof_lon[pos,-1], gprof_lat[pos,-1], color='red', lw = 0.1, )
    #plt.scatter(gprof_lon[pos,0], gprof_lat[pos,0], color='red', lw = 0.1, )


    plt.scatter(dy,dx, color='green', label='Deutschland Eckpunkte')

    plt.title(str(sc_y[pos]) +'-'+ str(sc_M[pos])+'-'+
              str(sc_dom[pos])+' '+str(sc_h[pos])+':'+
              str(sc_m[pos])+':'+str(sc_sec[pos])+ ' UTC')

    plt.legend(loc='lower left')
    plt.xlim(-180,180); plt.ylim(-90, 90)
    plt.xlabel('Lon'); plt.ylabel('Lat')
    plt.grid()
plt.show()

for jj in range(10):
    pos1 =  jj
    print str(sc_y[pos1]) +'-'+ str(sc_M[pos1])+'-'+str(sc_dom[pos1])+' '+str(sc_h[pos1])+':'+str(sc_m[pos1])+':'+\
          str(sc_sec[pos1])+ ' UTC,'+ ' H: ' + str(SC_alt[pos1]), 'Lon '+str(SC_lon[pos1])+ ' Orientation: '+ str(SC_ori[pos1])


'''
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
ax.plot(SC_lon, SC_lat, SC_alt)
plt.show()


#Given longitude and latitude on a sphere of radius R,
#the 3D coordinates P = (P.x, P.y, P.z) are:
#  P.x = R * cos(latitude) * cos(longitude)
#  P.y = R * cos(latitude) * sin(longitude)
#  P.z = R * sin(latitude)


radius = 6371.0
ri = radius + SC_alt
px,py,pz = ri * np.cos(SC_lat) *np.cos(SC_lon),\
           ri * np.cos(SC_lat) *np.sin(SC_lon),\
           ri * np.sin(SC_lat)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = Axes3D(fig)
ax.plot(py, px, pz, '-k')
plt.show()

'''
#Todo: Index Where is year == , m == d ==

#np.where(sc_y==2014)
#np.where(sc_m==30)
#np.where(sc_dom==21)
#np.where(sc_h==6)

