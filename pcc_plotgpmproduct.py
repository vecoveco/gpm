"""

Einlesen und darstellen von GPM Dateien


"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
import wradlib as wrl
from osgeo import osr
import os



ZP = '20141007023500'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'#'20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]


## Read GPROF
pfad2 = ('/home/velibor/shkgpm/data/example/gprof/*.HDF5')

pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]
gpmdprs = h5py.File(pfad_gprof_g, 'r')

gprof_lat=np.array(gpmdprs['S1']['Latitude'])			#(7934, 24)
gprof_lon=np.array(gpmdprs['S1']['Longitude'])			#(7934, 24)
gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])
gprof_pp[gprof_pp<0]=np.nan


## Read DPR
pfad2 = ('/home/velibor/shkgpm/data/example/dpr/*.HDF5')

pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]
gpmdprs = h5py.File(pfad_gprof_g, 'r')

dpr_lat=np.array(gpmdprs['NS']['Latitude'])			#(7934, 24)
dpr_lon=np.array(gpmdprs['NS']['Longitude'])			#(7934, 24)
dpr_pp=np.array(gpmdprs['NS']['SLV']['zFactorCorrectedNearSurface'])
dpr_pp[dpr_pp<0]=np.nan


plt.pcolormesh(gprof_lon, gprof_lat, np.ma.masked_invalid(gprof_pp), cmap='jet')
plt.plot(gprof_lon[:,0], gprof_lat[:,0], color='black')
plt.plot(gprof_lon[:,-1], gprof_lat[:,-1], color='black')
plt.plot(gprof_lon[:,87], gprof_lat[:,87], color='black')
plt.plot(gprof_lon[:,135], gprof_lat[:,135], color='black')

#plt.pcolormesh(dpr_lon, dpr_lat, np.ma.masked_invalid(dpr_pp), alpha=0.1, cmap='Blues')
plt.plot(dpr_lon[:,0], dpr_lat[:,0], color='red')
plt.plot(dpr_lon[:,-1], dpr_lat[:,-1], color='red')
plt.show()