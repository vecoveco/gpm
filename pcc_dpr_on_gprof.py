
import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats, linspace
import wradlib as wrl
from osgeo import osr
from pcc import get_time_of_gpm
from pcc import cut_the_swath

## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders
from pcc import plot_radar
from pcc import get_miub_cmap
my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()



ZP = '20141007'
year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

ye = ZP[2:4]

## Read GPM Data
## -------------

pfad = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gpmd = glob.glob(pfad)
pfad_gpmd_g = pfad_gpmd[0]

gpmdprd = h5py.File(pfad_gpmd_g, 'r')
dpr_lat = np.array(gpmdprd['NS']['Latitude'])
dpr_lon = np.array(gpmdprd['NS']['Longitude'])
#dpr_pp = np.array(gpmdprd['NS']['SLV']['zFactorCorrectedNearSurface'])
dpr_pp = np.array(gpmdprd['NS']['SLV']['precipRateNearSurface'])

dpr_pp[dpr_pp==-9999.9]= np.nan


pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
print pfad2
pfad_gprof = glob.glob(pfad2)
print pfad_gprof
pfad_gprof_g = pfad_gprof[0]

gpmdprs = h5py.File(pfad_gprof_g, 'r')
gprof_lat=np.array(gpmdprs['S1']['Latitude'])
gprof_lon=np.array(gpmdprs['S1']['Longitude'])
gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])
gprof_pp[gprof_pp==-9999.9]= np.nan

## Cut the GPM Swath
## ------------------

from pcc import cut_the_swath
glon, glat, gpp = cut_the_swath(gprof_lon,gprof_lat,gprof_pp, eu=True)
dlon, dlat, dpp = cut_the_swath(dpr_lon,dpr_lat,dpr_pp, eu=True)


proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

g_x, g_y = wradlib.georef.reproject(glon, glat, projection_target=proj_stereo , projection_source=proj_wgs)
g_xy = np.vstack((g_x.ravel(), g_y.ravel())).transpose()

d_x, d_y = wradlib.georef.reproject(dlon, dlat, projection_target=proj_stereo , projection_source=proj_wgs)
d_xy = np.vstack((d_x.ravel(), d_y.ravel())).transpose()



## INTERLOLATION
## --------------

gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gprof_xy = np.vstack((g_x.ravel(), g_y.ravel())).transpose()

grid_dpr_xy = np.vstack((d_x.ravel(), d_y.ravel())).transpose()


result = wrl.ipol.interpolate(grid_dpr_xy, grid_gprof_xy, dpp.reshape(dpp.shape[0]*dpp.shape[1],1), wrl.ipol.Idw, nnearest=4)

result = np.ma.masked_invalid(result)

rrr = result.reshape(g_x.shape)

plt.subplot(2,2,1)
plt.pcolormesh(g_x, g_y,np.ma.masked_invalid(gpp), vmin=0, vmax = 10)
plt.subplot(2,2,2)
plt.pcolormesh(d_x, d_y,np.ma.masked_invalid(dpp), vmin=0, vmax=10)
plt.subplot(2,2,3)
plt.pcolormesh(g_x, g_y,np.ma.masked_invalid(rrr), vmin=0, vmax=10)
plt.plot(d_x[:,0],d_y[:,0], color='black',lw=1)
plt.plot(d_x[:,-1],d_y[:,-1], color='black',lw=1)
plt.subplot(2,2,4)
rrr[rrr<=0]=np.nan
gpp[gpp<=0]=np.nan
plt.scatter(rrr, gpp)

plt.show()