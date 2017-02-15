"""

Dieses Programm soll die Interpolationsmethoden aus wradlib darstellen und verdeutlichen

"""

"""

Einlesen und darstellen von GPM DPR und Radolan Dateien

Radolanpfad:

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
#year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

ye = ZP[2:4]

## Read GPM Data
## -------------

pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gpm = glob.glob(pfad2)
pfad_gpm_g = pfad_gpm[0]

gpmdpr = h5py.File(pfad_gpm_g, 'r')
gprof_lat=np.array(gpmdpr['NS']['Latitude'])
gprof_lon=np.array(gpmdpr['NS']['Longitude'])

gprof_pp=np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])

gprof_pp[gprof_pp==-9999.9]= np.nan


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
#rwdata[rwdata < 0] = np.nan


## Cut the GPM Swath
## ------------------


blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


## INTERLOLATION
## --------------

gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

xy = np.vstack((x.ravel(), y.ravel())).transpose()

mask = ~np.isnan(rwdata)

rwdata_idw = rwdata.copy()
rwdata_near = rwdata.copy()
rwdata_lin = rwdata.copy()

# idw
result_idw = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata_idw.reshape(900*900,1),
                              wrl.ipol.Idw, nnearest=4)

result_idw = np.ma.masked_invalid(result_idw)

rrr_idw = result_idw.reshape(gpm_x.shape)
print 'idw fertig'

#near
result_near = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata_near.reshape(900*900,1),
                              wrl.ipol.Nearest)

result_near  = np.ma.masked_invalid(result_near )

rrr_near  = result_near .reshape(gpm_x.shape)
print 'near fertig'

#lin
result_lin = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata_lin.reshape(900*900,1),
                              wrl.ipol.Linear)

result_lin = np.ma.masked_invalid(result_lin)

rrr_lin = result_lin.reshape(gpm_x.shape)
print 'lin fertig'


## Interpolation of the binary Grid
## ------------------------------
res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1),
                               wrl.ipol.Idw, nnearest=4)
res_bin = res_bin.reshape(gpm_x.shape)

res_bin[res_bin!=0]= 1 # Randkorrektur

rand_y_unten = -4658.6447242655722
rand_y_oben = -3759.6447242655722
rand_x_rechts = 375.5378330781441


rrr_idw[np.where(gpm_y < rand_y_unten)] = np.nan
rrr_idw[np.where(gpm_y > rand_y_oben)] = np.nan
rrr_idw[np.where(gpm_x > rand_x_rechts)] = np.nan

rrr_lin[np.where(gpm_y < rand_y_unten)] = np.nan
rrr_lin[np.where(gpm_y > rand_y_oben)] = np.nan
rrr_lin[np.where(gpm_x > rand_x_rechts)] = np.nan

rrr_near[np.where(gpm_y < rand_y_unten)] = np.nan
rrr_near[np.where(gpm_y > rand_y_oben)] = np.nan
rrr_near[np.where(gpm_x > rand_x_rechts)] = np.nan

res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan
res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

ggg = gprof_pp_b * res_bin

## Nur Niederschlagsrelevante
rrr_idw[rrr_idw<5]=np.nan
rrr_near[rrr_near<5]=np.nan
rrr_lin[rrr_lin<5]=np.nan

ggg[ggg<0]=np.nan





ff = 15
cc = 0.5
fig = plt.figure(figsize=(16,16))


ax2 = fig.add_subplot(331, aspect='equal')#------------------------------------

pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr_lin),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Interpoliert Linear: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax2)
plot_radar(bonnlon, bonnlat, ax2, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax3 = fig.add_subplot(332, aspect='equal')#------------------------------------

pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr_idw),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Interpoliert IDW: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax3)
plot_radar(bonnlon, bonnlat, ax3, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

ax4 = fig.add_subplot(333, aspect='equal')#------------------------------------

pm4 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr_near),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Interpoliert Nearest: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax4)
plot_radar(bonnlon, bonnlat, ax4, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

from pcc import plot_scatter
ax1 = fig.add_subplot(334, aspect='equal')
plot_scatter(ggg, rrr_lin)
plt.xlabel('dpr')
plt.ylabel('lin')
ax4 = fig.add_subplot(335, aspect='equal')
plot_scatter(ggg, rrr_idw)
plt.xlabel('dpr')
plt.ylabel('idw')
ax5 = fig.add_subplot(336, aspect='equal')
plot_scatter(ggg, rrr_near)
plt.xlabel('dpr')
plt.ylabel('near')

ax7 = fig.add_subplot(338, aspect='equal')
pm7 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ggg - rrr_idw),
                     cmap=my_cmap, vmin=-3, vmax=3,zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Diff DPR-IDW: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax7)
plot_radar(bonnlon, bonnlat, ax7, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

ax8 = fig.add_subplot(337, aspect='equal')
pm8 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ggg - rrr_lin),
                     cmap=my_cmap, vmin=-3, vmax=3, zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Diff DPR - LIN: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax8)
plot_radar(bonnlon, bonnlat, ax8, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)

ax9 = fig.add_subplot(339, aspect='equal')
pm9 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ggg - rrr_near),
                     cmap=my_cmap, vmin=-3, vmax=3,
                     zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Reflectivity [dBZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Diff DPR-NEAR: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax9)
plot_radar(bonnlon, bonnlat, ax9, reproject=True)
plt.grid(color='r')
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)



plt.tight_layout()
plt.savefig('/home/velibor/shkgpm/plot/gpm_ipol_radolan_'+ZP + '.png' )
plt.close()
#plt.show()

from pcc import test
print test(rrr_near)
print test(rrr_lin)
print test(rrr_idw)






