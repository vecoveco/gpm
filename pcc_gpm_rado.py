"""

Einlesen und darstellen von GPM und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats
import wradlib as wrl
from osgeo import osr



ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]


# Zeitstempel nach YYYYMMDDhhmmss
#ZP = '20141007023500'#'20161024232500'#'20150427223500' #'20161024232500'##'20160917102000'#'20160917102000'#'20160805054500'
#ZP = '20170203005500'
#ZP = '20141007023500'# '20140629150000'
#ZP = '20140629150000'
#ZP = '20140826221000'
#ZP = '20140921071000'
#ZP = '20161024232500'
#ZP = '20161001060000'
ZP = '20140609132500'

year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------

pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

print 'RADOLAN PFAD: ', pfad

try:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad)
except EnvironmentError:
    rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

#Binaere Grid
rn = rwdata.copy()
rn[rn != -9999] = 1
rn[rn == -9999] = 0


radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5


## Read GPROF
## ------------
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

gpmdprs = h5py.File(pfad_gprof_g, 'r')
gprof_lat=np.array(gpmdprs['NS']['Latitude'])
gprof_lon=np.array(gpmdprs['NS']['Longitude'])

gprof_pp=np.array(gpmdprs['NS']['SLV']['zFactorCorrectedNearSurface'])

gprof_pp[gprof_pp==-9999.9]= np.nan




from pcc import cut_the_swath
blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)


gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()



## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders

from pcc import get_miub_cmap
my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()

##################################################################INTERLOLATION
gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

xy = np.vstack((x.ravel(), y.ravel())).transpose()

mask = ~np.isnan(rwdata)


result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata[mask].reshape(900*900,1), wrl.ipol.Idw, nnearest=4)

result = np.ma.masked_invalid(result)

rrr = result.reshape(gpm_x.shape)


# Interpolation des Binaear Grids
res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1), wrl.ipol.Idw, nnearest=4)
res_bin = res_bin.reshape(gpm_x.shape)
res_bin[res_bin!=0]= 1 #Randkorrektur


rand_y_unten = -4658.6447242655722
rand_y_oben = -3759.6447242655722
rand_x_rechts = 375.5378330781441


rrr[np.where(gpm_y < rand_y_unten)] = np.nan
rrr[np.where(gpm_y > rand_y_oben)] = np.nan
rrr[np.where(gpm_x > rand_x_rechts)] = np.nan

res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan

#Z = wradlib.trafo.idecibel(rwdata)
#rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


#Zr = wradlib.trafo.idecibel(rrr)
#rrr = wradlib.zr.z2r(Zr, a=200., b=1.6)
#
#
#
#
#rrr[rrr<=0]=0




from pcc import plot_radar

########################################################################## PLOT
###########################################################################----

ff = 15
cc = 0.5
fig = plt.figure(figsize=(14,10))
plt.suptitle('Problem: Changing RADOLAN observation section', fontsize=ff)
ax1 = fig.add_subplot(131, aspect='equal')
plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm1 = plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)

plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=cc)
cb.set_label("Ref [dbZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)



plot_borders(ax1)

plot_radar(bonnlon, bonnlat, ax1, reproject=True)

plt.title('RADOLAN Ref: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
#plt.xlim(-1000, 850)
#plt.ylim(-5500, -3000)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')

ax2 = fig.add_subplot(132, aspect='equal')
plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Ref [dbZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR Ref: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff)
plot_borders(ax2)
plot_radar(bonnlon, bonnlat, ax2, reproject=True)

plt.grid(color='r')
plt.tight_layout()

#plt.xlim(-1000, 850)
#plt.ylim(-5500, -3000)
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')


ax2 = fig.add_subplot(133, aspect='equal')
plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Ref [dbZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plt.title('RADOLAN Ref Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax2)
plot_radar(bonnlon, bonnlat, ax2, reproject=True)

plt.grid(color='r')
plt.tight_layout()

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')


plt.tight_layout()
plt.show()




##############################################################################

#Test Fig
ff = 15
fig = plt.figure(figsize=(10,10))

ax21 = fig.add_subplot(221, aspect='equal')

plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm21 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Ref Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax21)
plot_radar(bonnlon, bonnlat, ax21, reproject=True)

plt.grid(color='r')


ax44 = fig.add_subplot(222, aspect='equal')
plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm44 = plt.pcolormesh(gpm_x, gpm_y,res_bin,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2, alpha=0.2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Binaeres Grid: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax44)
plot_radar(bonnlon, bonnlat, ax44, reproject=True)

#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
plt.grid(color='r')


ax22 = fig.add_subplot(223, aspect='equal')
pm31 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR without BinGrid : \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax22)
plot_radar(bonnlon, bonnlat, ax22, reproject=True)

#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
plt.grid(color='r')


ax444 = fig.add_subplot(224, aspect='equal')
pm444 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b*res_bin),
                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Ref [dbz]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('GPM DPR with BinGrid: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax444)
plot_radar(bonnlon, bonnlat, ax444, reproject=True)

#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.show()



####Testplot2

ff = 15
fig = plt.figure(figsize=(10,10))

ax212 = fig.add_subplot(121, aspect='equal')

plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
#pm212 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
#                     cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)

#cb = plt.colorbar(shrink=0.8)
#cb.set_label("Ref [dbz]",fontsize=ff)
#cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Binaeres Grid',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax212)
plot_radar(bonnlon, bonnlat, ax212, reproject=True)

plt.grid(color='r')


ax442 = fig.add_subplot(122, aspect='equal')
plt.pcolormesh(x, y,rn*0.2,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2)
pm442 = plt.pcolormesh(gpm_x, gpm_y,res_bin,
                     cmap=my_cmap, vmin=0, vmax=1,zorder=2, alpha=0.2)

#cb = plt.colorbar(shrink=0.8)
#cb.set_label("Ref [dbz]",fontsize=ff)
#cb.ax.tick_params(labelsize=ff)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
#plt.xlabel("x [km] ",fontsize=ff)
#plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Binaeres Grid Interpolated',fontsize=ff) #RW Product Polar Stereo
plot_borders(ax442)
plot_radar(bonnlon, bonnlat, ax442, reproject=True)

#plt.xlim(-420,390)
#plt.ylim(-4700, -3700)
plt.grid(color='r')
plt.show()



