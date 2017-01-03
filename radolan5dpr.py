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
#import mpl_toolkits.basemap.pyproj as pyproj


ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
TH_rain= 0.2

# Zeitstempel nach YYYYMMDDhhmmss
ZP = '20160917102000'#'20161024232500'#'20140609132500'#'20160917102000'#'20160917102000'#'20160805054500'#'20141007023500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------
iii = 0
#pfad = ('/user/velibor/SHKGPM/data/radolan/*bin')
pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')
#print pfad
#pfad_radolan= sorted(glob.glob(pfad))
#print pfad_radolan
#pfad_radolan = pfad_radolan[iii]
pfad_radolan = pfad[:-3]

####### pfad

rw_filename = wradlib.util.get_wradlib_data_file(pfad)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)


rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

#sec = rwattrs['secondary']
#rwdata.flat[sec] = -9999
#rwdata = np.ma.masked_equal(rwdata, -9999)
radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]
#Z = wradlib.trafo.idecibel(rwdata)
#rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)

## Radolan lat lon Koordinaten
# ------------------------------
'''
## mask data
sec = rwattrs['secondary']
rwdata.flat[sec] = -9999
rwdata = np.ma.masked_equal(rwdata, -9999)

## create radolan projection object
proj_stereo = wradlib.georef.create_osr("dwd-radolan")

## create wgs84 projection object
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

## get radolan grid
radolan_grid_xy = wradlib.georef.get_radolan_grid(900, 900)
x1 = radolan_grid_xy[:, :, 0]
y1 = radolan_grid_xy[:, :, 1]

## convert to lonlat
radolan_grid_ll = wradlib.georef.reproject(radolan_grid_xy,
                                       projection_source=proj_stereo,
                                       projection_target=proj_wgs)
lon1 = radolan_grid_ll[:, :, 0]
lat1 = radolan_grid_ll[:, :, 1]
'''

#chose dpr corra gprof
ip = 0

GPMI = ['corra', 'dpr', 'gprof']
GPMI = GPMI[ip]
GPMI_name = ['CORRA RR (mm/h)', 'DPR RR (mm/h)','GPROF RR (mm/h)']

## Read GPROF
## ------------
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/'+str(GPMI)+'/*.HDF5')
pfad_gprof = glob.glob(pfad2)
print pfad_gprof
pfad_gprof_g = pfad_gprof[0]


gpmdprs = h5py.File(pfad_gprof_g, 'r')

if ip==0:
    gprof_lat = np.array(gpmdprs['NS']['Latitude'])
    gprof_lon = np.array(gpmdprs['NS']['Longitude'])
    gprof_pp = np.array(gpmdprs['NS']['surfPrecipTotRate'])

elif ip == 1:
    gprof_lat = np.array(gpmdprs['NS']['Latitude'])
    gprof_lon = np.array(gpmdprs['NS']['Longitude'])
    gprof_pp = np.array(gpmdprs['NS']['SLV']['precipRateNearSurface'])

elif ip == 2:
    gpmgmi_S1 = gpmdprs['S1']
    gprof_lat = np.array(gpmgmi_S1['Latitude'])
    gprof_lon = np.array(gpmgmi_S1['Longitude'])
    gprof_pp = np.array(gpmgmi_S1['surfacePrecipitation'])
    gprof_pp[gprof_pp<=0] = np.nan
else:
    'ip Falsch!'



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


alon = gprof_lon[latstart:latend]
alat = gprof_lat[latstart:latend]
gprof_pp_a = gprof_pp[latstart:latend]


ailat= np.where((alat>bonn_lat1) & (alat<bonn_lat2))
ailon= np.where((alon>bonn_lon1) & (alon<bonn_lon2))
alonstart = ailon[0][0]
alonend = ailon[0][-1]
alatstart = ailat[0][0]
alatend = ailat[0][-1]

blon = alon[alonstart:alonend]
blat = alat[alonstart:alonend]
gprof_pp_b = gprof_pp_a[alonstart:alonend]

print 'gprof min max:' + str(np.nanmin(gprof_pp_b)), str(np.nanmax(gprof_pp_b))
## GPM lon/lat in GK
## ------------
#proj_gk = osr.SpatialReference()
#proj_gk.ImportFromEPSG(31466)
#proj_ll = osr.SpatialReference()
#proj_ll.ImportFromEPSG(4326)
#gk3 = wradlib.georef.epsg_to_osr(31467)

proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

#, ,np.ma.masked_invalid(gprof_pp[latstart:latend]
#gpm_x, gpm_y = wrl.georef.reproject(gprof_lon[latstart:latend], gprof_lat[latstart:latend], projection_source=proj_ll,projection_target=proj_gk)
gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()





## Landgrenzenfunktion
## -------------------
'''
def plot_ocean(ax):
    # open the input data source and get the layer
    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    ocean, keys = wradlib.georef.get_shape_coordinates(inLayer)
    wradlib.vis.add_lines(ax, ocean, color='black', lw=2 , zorder=4)
'''
def plot_ocean(ax):

    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    inLayer.SetSpatialFilterRect(1, 45, 19, 56.5)
    borders, keys = wradlib.georef.get_shape_coordinates(inLayer)
    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)
    print borders.shape

    for j in range(borders.shape[0]):
        bu = np.array(borders[j].shape)

        a = np.array(bu.shape)

        if a==1:
            for i in range(0,borders[j].shape[0],1):
                bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
        if a==2:    #richtig
            bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

        bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

        wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
    ax.autoscale()

def plot_borders(ax):

    from osgeo import osr
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    india = osr.SpatialReference()
    # asia south albers equal area conic
    india.ImportFromEPSG(102028)

    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    # country list
    countries = ['Germany']#,'France','Denmark', 'Netherlands', 'Poland']
    # open the input data source and get the layer
    filename = wradlib.util.get_wradlib_data_file('/automount/db01/python/data/NED/10m/cultural/10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    # iterate over countries, filter accordingly, get coordinates and plot
    for item in countries:
        #print item
        # SQL-like selection syntax
        fattr = "(name='"+item+"')"
        inLayer.SetAttributeFilter(fattr)
        # get borders and names
        borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')

        for j in range(borders.shape[0]):
            bu = np.array(borders[j].shape)
            a = np.array(bu.shape)

            if a==1:
                for i in range(0,borders[j].shape[0],1):
                    bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                    bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                    wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
            if a==2:    #richtig
                bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
    ax.autoscale()
'''
def plot_dem(ax):
    filename = wradlib.util.get_wradlib_data_file('/home/velibor/cosmo_de_4326.tif')
    # pixel_spacing is in output units (lonlat)
    rastercoords, rastervalues = wradlib.io.read_raster_data(filename,
                                                         spacing=0.005)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                        rastervalues, cmap=plt.cm.gist_earth, norm=LogNorm(),
                        vmin=1, vmax=3000, zorder=1)
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    #cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    #cb = plt.gcf().colorbar(dem, cax=cax1,
                           #ticks=ticker.LogLocator(subs=range(10)))
    #cb.set_label('terrain height [m]')

def get_miub_cmap ():
    startcolor = 'white' # a dark olive
    color1 = '#8ec7ff' #'cyan' # a bright yellow
    color2 = 'dodgerblue'
    color3 = 'lime'
    color4 = 'yellow'
    color5 = 'darkorange'
    color6 = 'red'
    color7 = 'purple'
    #color6 = 'grey'
    endcolor = 'darkmagenta' # medium dark red
    colors = [startcolor, color1, color2, color3, color4, color5, color6, endcolor]
    return col.LinearSegmentedColormap.from_list('miub1' ,colors)
'''

dataset1, inLayer1 = wradlib.io.open_shape('/automount/db01/python/data/ADM/germany/vg250_0101.gk3.shape.ebenen/vg250_ebenen/vg250_l.shp')

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')


##################################################################INTERLOLATION
gk3 = wradlib.georef.epsg_to_osr(31467)

grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose() # GPM Grid erschaffen

xy = np.vstack((x.ravel(), y.ravel())).transpose()

mask = ~np.isnan(rwdata)

result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata[mask].reshape(900*900,1), wrl.ipol.Idw, nnearest=4)  #Idw

result = np.ma.masked_invalid(result)

rrr = result.reshape(gpm_x.shape)

#Todo:  FEHLER IWO...wahrscheinlich Randbedingungen!? Fehler Wert 10
#SCHEMA http://wradlib.org/wradlib-docs/latest/notebooks/interpolation/wradlib_ipol_example.html

Z = wradlib.trafo.idecibel(rwdata)
rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)


Zr = wradlib.trafo.idecibel(rrr)
rrr = wradlib.zr.z2r(Zr, a=200., b=1.6)



print 'rwdata min max:' + str(np.nanmin(rwdata)), str(np.nanmax(rwdata))

print 'rrr min max:' + str(np.nanmin(rrr)), str(np.nanmax(rrr))

## PLOT
## ----
ff = 15
mini, maxi = 0.1, 10
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(221, aspect='equal')
plt.pcolormesh(x, y, rwdata, cmap=my_cmap,vmin=mini,vmax=maxi, zorder=2)
#plt.scatter(x, y, rwdata, cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plot_borders(ax1)
plt.title('RADOLAN Rainrate: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo

#plot_ocean(ax1)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)
plt.grid(color='r')
plt.xlim(-420,390)
plt.ylim(-4700, -3700)


ax2 = fig.add_subplot(224, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                     cmap=my_cmap,vmin=mini,vmax=maxi, zorder=2)

#pm2 = plt.pcolormesh(gprof_lon[latstart:latend], gprof_lat[latstart:latend],np.ma.masked_invalid(gprof_pp[latstart:latend]),
                     #cmap=my_cmap,vmin=0.1,vmax=10, zorder=2)
cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title(GPMI_name[ip] + ': \n' + str(pfad_gprof_g[66:70]) + '-' +str(pfad_gprof_g[70:72])+ '-' +
          str(pfad_gprof_g[72:74]) + ' T: ' +str(pfad_gprof_g[76:82]) + '-' + str(pfad_gprof_g[84:90]) + ' UTC',fontsize=ff)
#plot_ocean(ax2)
plot_borders(ax2)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(ax1.get_ylim())
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)




ax2 = fig.add_subplot(223, aspect='equal')
pm2 = plt.pcolormesh(gpm_x, gpm_y,rrr,
                     cmap=my_cmap,vmin=mini,vmax=maxi, zorder=2)

cb = plt.colorbar(shrink=0.8)
cb.set_label("Rainrate (mm/h)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.xlabel("x [km] ",fontsize=ff)
plt.ylabel("y [km]  ",fontsize=ff)
plt.title('RADOLAN Rainrate Interpoliert: \n'+'20' + str(pfad_radolan[-20:-18])+'-'+str(pfad_radolan[-18:-16])+'-'+str(pfad_radolan[-16:-14])+
       ' T: '+str(pfad_radolan[-14:-10]) + '00 UTC',fontsize=ff) #RW Product Polar Stereo
#plot_ocean(ax2)
plot_borders(ax2)
#plt.xticks(fontsize=ff)
#plt.yticks(fontsize=ff)
#plt.xlim((bonn_lon1-1,bonn_lon2+1))
#plt.ylim((bonn_lat1-1,bonn_lat2+1))
plt.grid(color='r')
plt.tight_layout()
#Limits Setzen
ax2.set_xlim(ax1.get_xlim())
ax2.set_ylim(ax1.get_ylim())
#plt.xticks(fontsize=0)
#plt.yticks(fontsize=0)


ax2 = fig.add_subplot(222, aspect='equal')

A = rrr
B = np.ma.masked_invalid(gprof_pp_b)

ref = np.copy(rrr)
est = np.copy(np.ma.masked_invalid(gprof_pp_b))

A[A<TH_rain] = np.nan
B[B<TH_rain] = np.nan



mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept

##############################################################################
#plt.hist2d(B[mask],A[mask], bins=100, norm=LogNorm())
#plt.hexbin(B[mask],A[mask])
#plt.colorbar()
#plt.show()

xx = B[mask]
yy = A[mask]
xedges, yedges = np.linspace(-4, 4, 42), np.linspace(-25, 25, 42)
hist, xedges, yedges = np.histogram2d(xx, yy, (xedges, yedges))
xidx = np.clip(np.digitize(xx, xedges), 0, hist.shape[0]-1)
yidx = np.clip(np.digitize(yy, yedges), 0, hist.shape[1]-1)
c = hist[xidx, yidx]
plt.scatter(xx, yy, c=c, label='RR [mm/h]')
plt.colorbar()
plt.plot(B,line,'r-')
maxAB = np.nanmax([np.nanmax(xx),np.nanmax(yy)])
#log or linear
plt.yscale('linear')
plt.xscale('linear')
plt.xlim(0,maxAB + 1)
plt.ylim(0,maxAB + 1)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=2, fancybox=True, shadow=True,
                    fontsize='small', title= "Slope: " + str(round(slope,3))
                                            + ', Intercept: '+  str(round(intercept,3)) + "\n Correlation: " +
                                            str(round(r_value,3)) + ', Std_err: '+  str(round(std_err,3)))
plt.xlabel(GPMI_name[ip])
plt.ylabel("RADOLAN RR [mm/h]")
plt.title(" .")

plt.grid(True)
plt.tight_layout()
plt.show()
###############################################################################

import pcc
R = pcc.skill_score(est,ref,0.1)
pcc.plot_score(est,ref,R)
plt.xlabel(GPMI_name[ip])
plt.ylabel("RADOLAN RR [mm/h]")
#plt.yscale('log')
#plt.xscale('log')
plt.title(" .")
plt.tight_layout()
plt.show()



'''
plt.hist(A[mask],bins=200, color='red', alpha=0.4, label='RADOLAN interpoliert')

plt.hist(B[mask], bins=200, color='blue', alpha=0.4, label='GPROF')
#pdf = stats.norm.pdf(sorted(B[mask]), B[mask], B[mask])
#plt.plot(sorted(B[mask]), pdf)
plt.xlabel("RR [mm/h]")
plt.ylabel("frequency")
plt.legend(loc='upper right')
plt.grid()
plt.show()


RR1_3 = A
GR1_3 = B
# forschleife in ein scatter mit verschiedenen farben
colors = ['blue', 'red','green']
rmin = [1,5,9]
rmax = [5,9,40]

for ii in range(3):
    Rmin, Rmax = rmin[ii], rmax[ii]

    RR1_3[RR1_3 < Rmin] = np.nan
    GR1_3[GR1_3 > Rmax] = np.nan

    RR1_3[RR1_3 < Rmin] = np.nan
    GR1_3[GR1_3 > Rmax] = np.nan

    mask = ~np.isnan(GR1_3) & ~np.isnan(RR1_3)
    slope, intercept, r_v, p_value, std_err = stats.linregress(GR1_3[mask], RR1_3[mask])
    line = slope*GR1_3+intercept

    plt.scatter(RR1_3,GR1_3, color=colors[ii], label=str(rmse(RR1_3,GR1_3))+'___'+str(r_v))
    plt.plot(GR1_3,line,color=colors[ii])


    from pcc import rmse

plt.legend(loc='lower right')
plt.grid()
plt.show()


'''