"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

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




ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
TH_rain= 0.2

# Zeitstempel nach YYYYMMDDhhmmss
ZP = '20160917102000'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'#'20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]



################################################################### Read GPROF
## ------------
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

gpmgmi = h5py.File(pfad_gprof_g, 'r')

gpmgmi.keys()
gpmgmi_S1=gpmgmi['S1']
gprof_lat=np.array(gpmgmi_S1['Latitude'])
gprof_lon=np.array(gpmgmi_S1['Longitude'])
gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
gprof_pp[gprof_pp<=0] = np.nan


bonn_lat1 = 47.9400
bonn_lat2 = 55.3500
bonn_lon1 = 6.40000
bonn_lon2 = 14.10000

ilat= np.where((gprof_lat>bonn_lat1) & (gprof_lat<bonn_lat2))
ilon= np.where((gprof_lon>bonn_lon1) & (gprof_lon<bonn_lon2))
lonstart = ilon[0][0]
lonend = ilon[0][-1]
latstart = ilat[0][0]
latend = ilat[0][-1]


############################################################ GPM lon/lat in GK
## ------------
proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(gprof_lon[latstart:latend],
                                        gprof_lat[latstart:latend],
                                        projection_target=proj_stereo ,
                                        projection_source=proj_wgs)

grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

