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



# Zeitstempel nach YYYYMMDDhhmmss
#'20161024232500'# #'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'
#ZP = '20170203005500'
#ZP = '20141007023500'# '20140629150000'
#ZP = '20140629150000'
#ZP = '20140826221000'
ZP = '20140921071000'
#ZP = '20161024232500'
#ZP = '20161001060000'
#ZP = '20150427223500'

year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]


## Read RADOLAN Data
## -----------------

r_pro = 'rz'

pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
        str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
        str(ye)+str(m)+ str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

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
rwdata = np.ma.masked_equal(rwdata, -9999) *8
#rwdata[rwdata < 0] = np.nan



## Read GPM Data
## -------------

pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad2)
pfad_gprof_g = pfad_gprof[0]

gpmdprs = h5py.File(pfad_gprof_g, 'r')
gprof_lat=np.array(gpmdprs['S1']['Latitude'])
gprof_lon=np.array(gpmdprs['S1']['Longitude'])

gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])

gprof_pp[gprof_pp==-9999.9]= np.nan



## Cut the GPM Swath
## ------------------

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
from pcc import plot_radar
from pcc import get_miub_cmap
my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()


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


rrr[np.where(gpm_y < rand_y_unten)] = np.nan
rrr[np.where(gpm_y > rand_y_oben)] = np.nan
rrr[np.where(gpm_x > rand_x_rechts)] = np.nan

res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan
res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

ggg = gprof_pp_b * res_bin

## Nur Niederschlagsrelevante
rrr[rrr<0]=np.nan
ggg[ggg<0]=np.nan





ff = 15
cc = 0.5
maxi = np.nanmax([ggg,rrr])
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(221, aspect='equal')

pm1 = plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)

plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Ref [dbZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)

plot_borders(ax1)
plot_radar(bonnlon, bonnlat, ax1, reproject=True)

plt.title('RADOLAN Ref: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff)
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

ax2 = fig.add_subplot(222, aspect='equal')

pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ggg),
                     cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
cb = plt.colorbar(shrink=cc)
cb.set_label("Ref [dbZ]",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.title('GPM DPR Ref: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff)
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

ax2 = fig.add_subplot(223, aspect='equal')

pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                     cmap=my_cmap, vmin=0.01, vmax=10,zorder=2)
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

ax4 = fig.add_subplot(224, aspect='equal')

maske = ~np.isnan(ggg) & ~np.isnan(rrr)
slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
line = slope * ggg +intercept

from pcc import skill_score
SS = skill_score(ggg,rrr,th=0)

plt.scatter(ggg, rrr, label='Reflectivity [dbZ]', color='grey')

text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
           '\ncorr: ' + str(round(r_value,3)) + r'$\pm$: '+  str(round(std_err,3))+
        '\nN: '+ str(SS['N'])+
        '\nHit: ' + str(round(SS['H']/SS['N'],3)*100)+'%'+
        '\nMiss: ' + str(round(SS['M']/SS['N'],3)*100)+'%'+
        '\nFalse: ' + str(round(SS['F']/SS['N'],3)*100)+'%'+
        '\nCnegative: ' + str(round(SS['C']/SS['N'],3)*100)+'%'+
        '\nPOD: ' + str(round(SS['POD'],3))+
        '\nFAR: ' + str(round(SS['FAR'],3))+
        '\nBID: ' + str(round(SS['BID'],3))+
        '\nHSS: ' + str(round(SS['HSS'],3))+
        '\nHSS: w.i.p' +
        '\nHSS: w.i.p'
        )

ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=10,
                horizontalalignment='left', verticalalignment='top')

t1 = linspace(0,maxi,maxi)
plt.plot(t1,t1,'k-')
plt.plot(t1,t1 + 5,'k-.')
plt.plot(t1,t1 - 5,'k-.')
plt.plot(t1, t1*slope + intercept, 'r-', lw=3 )

plt.xlim(0,maxi)
plt.ylim(0,maxi)


plt.xlabel('GPM DPR Reflectivity [dbZ]',fontsize=ff)
plt.ylabel('RADOLAN Reflectivity [dbZ]',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid(color='r')


plt.tight_layout()
plt.show()



print 'f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3))


