"""

Interpolation Methode To New Grid

"""

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


ZP = '20141007'
year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
ye = ZP[2:4]

## Read GPM Data
## -------------

pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gpm = glob.glob(pfad2)
pfad_gpm_g = pfad_gpm[0]

gpmdpr = h5py.File(pfad_gpm_g, 'r')
gprof_lat = np.array(gpmdpr['NS']['Latitude'])
gprof_lon = np.array(gpmdpr['NS']['Longitude'])

gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])
#gprof_pp[gprof_pp==-9999.9]= np.nan

gpm_time = gpmdpr['NS']['ScanTime']
gpm_zeit = get_time_of_gpm(gprof_lon, gprof_lat, gpm_time)

ht, mt = gpm_zeit[14:16], str(int(round(float(gpm_zeit[17:19])/5.0)*5.0))
if mt == '0':
    mt = '00'
elif mt == '5':
    mt = '05'
print mt
print gpm_zeit


## Read RADOLAN Data
## -----------------
r_pro = 'rx'

pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
        str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
        str(ye)+str(m)+ str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

try:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad)
except EnvironmentError:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
#Binaere Grid
rn = rwdata.copy()
rn[rn != -9999] = 1
rn[rn == -9999] = 0

radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5


## Cut the GPM Swath
## ------------------
blon, blat, bpp = cut_the_swath(gprof_lon,gprof_lat,gprof_pp, eu=0)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

gn = bpp.copy()
gn[gn != -9999.9] = 1
gn[gn == -9999.9] = 0

## INTERLOLATION
## --------------

gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

xy = np.vstack((x.ravel(), y.ravel())).transpose()

mask = ~np.isnan(rwdata)

result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=4)

result = np.ma.masked_invalid(result)

rrr = result.reshape(gpm_x.shape)

## Interpolation of the binary Grid
## ------------------------------
res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1), wrl.ipol.Idw, nnearest=4)
res_bin = res_bin.reshape(gpm_x.shape)

res_bin[res_bin!=0]= 1 # Randkorrektur

rand_y_unten = -4658.6447242655722
rand_y_oben = -3759.6447242655722
rand_x_rechts = 375.5378330781441


'''
rrr[np.where(gpm_y < rand_y_unten)] = np.nan
rrr[np.where(gpm_y > rand_y_oben)] = np.nan
rrr[np.where(gpm_x > rand_x_rechts)] = np.nan

res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan
res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

ggg = bpp * res_bin

## Dynamischer Threshold
THref = np.nanmax([np.nanmin(rrr),np.nanmin(ggg)])

## Nur Niederschlagsrelevante
rrr[rrr < THref]=np.nan
ggg[ggg < THref]=np.nan

gn = ggg.copy()
where_are_NaNs = np.isnan(gn)
gn[where_are_NaNs] = 0
gn[gn>0]=1

'''
## INTERLOLATION to new grid
## --------------------------
grid_size = 10

x_new, y_new = np.arange(-500,400,grid_size), np.arange(-4660, -3760,grid_size)

xx_new, yy_new = np.meshgrid(x_new, y_new)

xy_new = np.vstack((xx_new.ravel(), yy_new.ravel())).transpose()


bpp = np.ma.masked_equal(bpp, -9999)

result_gpm = wrl.ipol.interpolate(grid_gpm_xy, xy_new,  bpp.reshape(bpp.shape[0]*bpp.shape[1],1),
                                  wrl.ipol.Idw, nnearest=4)
result = np.ma.masked_invalid(result)
result_gpm = result_gpm.reshape(xx_new.shape)


result_rad = wrl.ipol.interpolate(xy, xy_new,  rwdata.reshape(rwdata.shape[0]*rwdata.shape[1],1),
                                  wrl.ipol.Idw, nnearest=4)
result_rad = result_rad.reshape(xx_new.shape)


res_gn = wrl.ipol.interpolate(grid_gpm_xy, xy_new, res_bin.reshape(res_bin.shape[0]*res_bin.shape[1],1),
                                  wrl.ipol.Idw, nnearest=4)
res_gn = res_gn.reshape(xx_new.shape)
res_gn=res_gn/res_gn


res_gpm = wrl.ipol.interpolate(grid_gpm_xy, xy_new, gn.reshape(gn.shape[0]*gn.shape[1],1),
                                  wrl.ipol.Idw, nnearest=4)
res_gpm = res_gpm.reshape(xx_new.shape)
res_gpm[res_gpm!=0]= 1 # Randkorrektur



from pcc import get_my_cmap

#result_gpm[np.where(xx_new < gpm_x[:,0])] = np.nan
#result_gpm[np.where(yy_new < gpm_y[:,-1])] = np.nan


plt.subplot(2,3,1)
plt.pcolormesh(gpm_x, gpm_y, np.ma.masked_invalid(rrr), cmap=get_my_cmap(), vmin=0.01, vmax=50)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)

plt.subplot(2,3,2)
plt.pcolormesh(gpm_x, gpm_y, np.ma.masked_invalid(bpp), cmap=get_my_cmap(), vmin=0.01, vmax=50)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
plt.subplot(2,3,4)
plt.pcolormesh(xx_new, yy_new, result_rad, cmap=get_my_cmap(), vmin=0.01, vmax=50)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)

plt.subplot(2,3,5)
plt.pcolormesh(xx_new, yy_new, np.ma.masked_invalid(result_gpm*res_gpm), cmap=get_my_cmap(), vmin=0.01, vmax=50)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)

plt.subplot(2,3,3)
A = result_gpm*res_gn
B = result_rad*res_gn

A[A<0]=np.nan
B[B<0]=np.nan



maske = ~np.isnan(A) & ~np.isnan(B)
slope, intercept, r_value, p_value, std_err = stats.linregress(A[maske], B[maske])
line = slope * A +intercept
plt.scatter(A,B)
plt.title(str(r_value))

plt.subplot(2,3,6)
plt.pcolormesh(xx_new, yy_new, res_gpm, cmap=get_my_cmap(), vmin=0.01, vmax=50)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)

plt.show()


