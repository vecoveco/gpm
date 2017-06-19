#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
# ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/

import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')
from pcc import get_miub_cmap as my_cmap
from pcc import plot_radar
from pcc import boxpol_pos
import wradlib as wrl
from osgeo import osr

Pos = boxpol_pos()
bblon, bblat = Pos['lon_ppi'], Pos['lat_ppi']
gkx, gky = Pos['gkx_ppi'], Pos['gky_ppi']

def get_frist_and_last_valid(data,attributes,thr):
    '''
    This function enables you to find valid pixels with preciptitation. To choose correct and uncluttered bins different
    radii are treatet differently:

        < 5 km all pixel must be larger than rhohv and zh threshold within a radius of +- 3 pixel
        >=5 km but < 20 km all pixel must be larger than rhohv and zh threshold within a radius of +- 2 pixel
        >=20 km but < 50 km all pixel must be larger than rhohv and zh threshold within a radius of +- 1 pixel
        >=50 km but < 100 km all pixel must be larger than rhohv and zh threshold within a bin_radius in equal azi of +- 1 pixel

    The first valid and the last valid radii are choosen as the ones where all substractions are fullfilled.
    '''

    zh_extend    = np.vstack(( data['SCAN0']['ZH']['data'][-3:,:], data['SCAN0']['ZH']['data'],
                              data['SCAN0']['ZH']['data'][:3,:]))
    zh_extend = zh_extend+2
    rhohv_extend = np.vstack(( data['SCAN0']['RHOHV']['data'][-3:,:], data['SCAN0']['RHOHV']['data'],
data['SCAN0']['RHOHV']['data'][:3,:]))

    a,b = zh_extend.shape
    bin_range=attributes['SCAN0']['bin_range']/1000 # in km

    r_valid        = np.zeros((a,b))
    r_start, r_end = np.zeros((a)), np.zeros((a))
    for j in range(3,a-3):
        for i in range(4,b):
            if   i*bin_range>50 and i*bin_range<100 and (zh_extend[j,i-1:i+2]>=thr['ZH']).all()       and (rhohv_extend[j,i-1:i+2]>=thr['rhohv']).all():
                r_valid[j,i]=1.
            elif i*bin_range>20 and i*bin_range<=50 and (zh_extend[j-1:j+2,i-1:i+2]>=thr['ZH']).all() and (rhohv_extend[j-1:j+2,i-1:i+2]>=thr['rhohv']).all():
                r_valid[j,i]=1.
            elif i*bin_range>5  and i*bin_range<=20 and (zh_extend[j-2:j+3,i-2:i+3]>=thr['ZH']).all() and (rhohv_extend[j-2:j+3,i-2:i+3]>=thr['rhohv']).all():
                r_valid[j,i]=1.
            elif i*bin_range<=5                     and (zh_extend[j-3:j+4,i-3:i+4]>=thr['ZH']).all() and (rhohv_extend[j-3:j+4,i-3:i+4]>=thr['rhohv']).all():
                r_valid[j,i]=1.

            #if data['SCAN0']['ZH']['data'][j,i]>=thr['ZH'] and data['SCAN0']['RHOHV']['data'][j,i]>=thr['rhohv']:
            #    r_valid[j,i]=1.
        if True in (i>0 for i in r_valid[j,:]):
            ri         = np.where(r_valid[j,:]>0)
            r_start[j] = ri[0][0]
            r_end[j]   = ri[0][-1]

    r_start = np.array(r_start[3:-3])
    r_end   = np.array(r_end[3:-3])

    return(r_start,r_end)

################################################################################################################

def smoothed_phidp(data, attributes, method='median'):

    '''
    This function enables you to smooth phidp via different methods (in theory). Yet there will only be one
    method implemented as median, maybe also mean. median will be default.
    '''

    a,b = data['SCAN0']['PHIDP']['data'].shape
    phidp = np.zeros((a,b))
    delta=10

    for j in range(a):
        for i in range(4,b):
            if i<=b-delta:
                phidp[j,i] = np.median(data['SCAN0']['PHIDP']['data'][j,i-delta:i+delta+1])
            else:
                phidp[j,i] = np.median(data['SCAN0']['PHIDP']['data'][j,i-delta:-1])

    return(phidp)

################################################################################################################

def ZPHI_Method(data,attributes,alpha=0.28):
    '''
    Calculations for ZPHI_Method:


    A(r)=[ z^b * C(b,PIA)] / [ I(r1,r2) + C(b,PIA) * I(r,r2) ]


    With:

    r1 < r < r2

    I(r1,r2) = 0.46 * b * Int_r1^r2( z^b ds)

    I(r,r2) = 0.46 * b * Int_r^r2( z^b ds)

    C(b,PIA) = 10^(0.1*b*alpha*deltaPhiDP)-1
    DeltaphiDP = Phidp(r2) - PhiDP(r1)

    r1 and r2 are picked by the function get_first_and_last_valid()
    '''
    #Thresholds:
    thr = {'rhohv' : 0.85 ,'ZH' : 5.}

    # Constants:
    b   = 0.78

    ZH    = data['SCAN0']['ZH']['data'] #shape(360,1000)
    ZH=ZH + 2.
    zh    = 10.**(ZH/10.)
    phidp = smoothed_phidp(data, attributes)
    resolution = attributes['SCAN0']['bin_range']/1000.
    #print(resolution)

    r1,r2 = get_frist_and_last_valid(data, attributes, thr)

    Ir1r2      = np.zeros(zh.shape[0])
    CbPIA      = np.zeros(zh.shape[0])
    deltaphidp = np.zeros(zh.shape[0])
    Irr2       = np.zeros(zh.shape)
    Ah         = np.empty(zh.shape)
    Ah.fill(np.nan)
    ZH_add     = np.empty(zh.shape)
    ZH_add.fill(np.nan)
    phidp_est  = np.empty(zh.shape)
    phidp_est.fill(np.nan)
    for j in range(zh.shape[0]):
        deltaphidp[j] = phidp[j,r2[j]]-phidp[j,r1[j]]
        if deltaphidp[j] >= 3.:
            Ir1r2[j]      = 0.46*b*np.sum(zh[j,r1[j]:r2[j]+1]**b)*resolution
            CbPIA[j]      = (10.**(0.1*alpha*b*deltaphidp[j]))-1.

            for i in np.arange(r1[j],r2[j]):
                Irr2[j,i] = 0.46*b*np.sum(zh[j,i:r2[j]+1]**b)*resolution

                Ah[j,i]        = ((zh[j,i]**b) * CbPIA[j]) / (Ir1r2[j] + (CbPIA[j] * Irr2[j,i]))
                Ah[Ah<0]       = np.nan
                ZH_add[j,i]    = 2*np.nansum(Ah[j,0:i+1])*resolution
                phidp_est[j,i] = (2./alpha)*np.nansum(Ah[j,0:i+1])*resolution
        else:
            Ah[j,:]        = np.nan
            ZH_add[j,:]    = np.nan
            phidp_est[j,:] = np.nan

        phidp_est[j,r2[j]:]=phidp_est[j,r2[j]]
    ZH_add[np.isnan(ZH_add)]=0.
    ZH_corr= ZH + ZH_add

    return (Ah,ZH,ZH_corr,phidp,phidp_est,deltaphidp,r1,r2)



##################################################################################################################
##################################################################################################################





'''Dieses Program soll dazu dienen die
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Reflektivitat zu validieren.
Hier werden mehrere Ueberflug analysiert'''



# Pfad mit String
# ---------------

# Hohe von DPR
TH = 12 #Threshold um Nullen fuer Niederschlag raus zu filtern

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]
offset = 2


#ZP = '20141007023744'
#ZP = '20140704134500'
#ZP = '20150225163500'
#ZP = '20150816070500'
#good#ZP = '20151208213500'
#good#ZP = '20151216024501'
#ZP = '20160503024500'
ZP = '20160209103500'

year = ZP[0:4]
m = ZP[4:6]
d = ZP[6:8]
ht = ZP[8:10]
mt = ZP[10:12]
st = ZP[12:14]



pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
print pfad_radar
#pfad_radar = sorted(glob.glob(pfad))
#print pfad_radar
pfad_radar_Ku = pfad_radar[0]

#/automount/radar-archiv/scans/2014/2014-10/2014-10-07/n_ppi_010deg

try:
    ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]
    #ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/n_ppi_280deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

except:
    ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",*.mvol")[0]

print ppi_datapath


# PPI BoxPol Daten einlesen
#---------------------------

ppi=h5py.File(ppi_datapath,'r')
data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

ZH = data['SCAN0']['ZH']['data']
PHIDP = data['SCAN0']['PHIDP']['data']
rho = data['SCAN0']['RHOHV']['data']

r = attrs['SCAN0']['r']
az = attrs['SCAN0']['az']
lon_ppi = attrs['VOL']['Longitude']
lat_ppi = attrs['VOL']['Latitude']
alt_ppi = attrs['VOL']['Height']


data,attributes=wrl.io.read_GAMIC_hdf5('/automount/radar-archiv/scans/2014/2014-10/2014-10-07/ppi_1p5deg/2014-10-07--02:37:44,00.mvol')
Ah,ZH,ZH_corr,phidp,phidp_est,deltaphidp,r1,r2 = ZPHI_Method(data,attributes)

R = ZH_corr


from wradlib.trafo import idecibel, decibel
R = idecibel(R)

#R = R + pia
R[151:165]=np.nan

#R[np.where(rho<0.95)]=np.nan
R[rho<=0.8]=np.nan

# DPR Einlesen
# ------------

gpmku = h5py.File(pfad_radar_Ku, 'r')
gpmku_HS=gpmku['NS']['SLV']
ku_lat=np.array(gpmku['NS']['Latitude'])			#(7934, 24)
ku_lon=np.array(gpmku['NS']['Longitude'])			#(7934, 24)
ku_pp=np.array(gpmku_HS['zFactorCorrectedNearSurface'])




# Lon Lat Bestimmung
# ------------------
radars = [ku_pp]
rad_lat = [ku_lat]
rad_lon = [ku_lon]
rad_name = ['DPR']

ii = 0

dpr_pp = radars[ii]
dpr_lat = rad_lat[ii]
dpr_lon = rad_lon[ii]
radarname = rad_name[ii]

bonn_lat1 = 49.9400
bonn_lat2 = 51.3500
bonn_lon1 = 6.40000
bonn_lon2 = 8.10000

ilat= np.where((dpr_lat>49.9400) & (dpr_lat<51.3500))
ilon= np.where((dpr_lon>6.40000) & (dpr_lon<8.10000))
lonstart = ilon[0][0]
lonend = ilon[0][-1]
latstart = ilat[0][0]
latend = ilat[0][-1]
dpr_pp[dpr_pp==-9999] = np.nan
#dpr_pp = dpr_pp[:,:,hi]  # Untersete Schicht

radar_location = (lon_ppi, lat_ppi, alt_ppi)
elevation = 1.5
azimuths = az
ranges = r
polargrid = np.meshgrid(ranges, azimuths)
lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

gk3 = wradlib.georef.epsg_to_osr(31467)
x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
xgrid, ygrid = wradlib.georef.reproject(dpr_lon[latstart:latend], dpr_lat[latstart:latend], projection_target=gk3)

grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()


xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[0],nnearest=40,p=2)
gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)

gridded = decibel(gridded)

R = decibel(R)

# ON RADOLAN GRID
proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)

box_x, box_y = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo , projection_source=proj_wgs)
gpm_x, gpm_y = wradlib.georef.reproject(dpr_lon[latstart:latend], dpr_lat[latstart:latend], projection_target=proj_stereo , projection_source=proj_wgs)

#Todo IM PLOT VERWENDEN !!!!!!!

# Plot
# ----



################################################################Swap!
#rrr, ggg = ggg, rrr

ff = 20
cc = 1.0
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------
#ax1, pm1 = wradlib.vis.plot_ppi(R,r,az,vmin=0.01,vmax=50, cmap=my_cmap())
pm1 = plt.pcolormesh(box_x, box_y, R, vmin=0, vmax=50, cmap=my_cmap())

#cb = plt.colorbar(pm1,shrink=cc)
#cb.set_label("Reflectivity (dBZ)",fontsize=ff)
#cb.ax.tick_params(labelsize=ff)

plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')

from pcc import plot_borders
plot_borders(ax1)
plot_radar(bblon, bblat, ax1, reproject=True, cband=False,col='black')

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')

#plt.title('BoXPol Reflectivity:\n 2014-10-07--02:37:44',fontsize=ff)

plt.grid(color='r')
plt.xlim(-335, -80)
plt.ylim(-4400,-4100)



ax2 = fig.add_subplot(222, aspect='equal')#------------------------------------

pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(dpr_pp[latstart:latend]),
                     cmap=my_cmap(), vmin=0.01, vmax=50, zorder=2)

plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')

plot_borders(ax2)
plot_radar(bblon, bblat, ax2, reproject=True, cband=False,col='black')


cb = plt.colorbar(pm2,shrink=cc)
cb.set_label("Reflectivity (dBZ)",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
#plt.title('GPM DPR Reflectivity \n 2014-10-07--02:35:27 ',fontsize=ff)
plt.grid(color='r')

plt.xlim(-335, -80)
plt.ylim(-4400,-4100)

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')


ax3 = fig.add_subplot(223, aspect='equal')#------------------------------------

pm3 = plt.pcolormesh(gpm_x, gpm_y, gridded,
                     cmap=my_cmap(), vmin=0.01, vmax=50,zorder=2)
plt.plot(gpm_x[:,0], gpm_y[:,0], color='black')
plt.plot(gpm_x[:,-1], gpm_y[:,-1], color='black')
plt.plot(gpm_x[:,23], gpm_y[:,23], color='black', ls='--')


#cb = plt.colorbar(pm3, shrink=cc)
#cb.set_label("Reflectivity (dBZ)",fontsize=ff)
#cb.ax.tick_params(labelsize=ff)

#plt.title('BoXPol Reflectivity Interpolated\n 2014-10-07--02:37:44',fontsize=ff) #RW Product Polar Stereo

plot_borders(ax3)
plot_radar(bblon, bblat, ax3, reproject=True, cband=False,col='black')

plt.grid(color='r')
plt.xlim(-335, -80)
plt.ylim(-4400,-4100)

plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')


ax4 = fig.add_subplot(224, aspect='equal')#------------------------------------

rrr = gridded.copy()
ggg = np.ma.masked_invalid(dpr_pp)[latstart:latend].copy()

rrr[rrr<TH]=np.nan
ggg[ggg<TH]=np.nan


maske = ~np.isnan(ggg) & ~np.isnan(rrr)
slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
line = slope * ggg +intercept

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(rrr[maske], ggg[maske])
line2 = slope2 * rrr +intercept2

diffi = ggg[maske]-rrr[maske]
bias = np.nansum(diffi)/len(diffi)
rmse = np.sqrt(np.nansum(((diffi)**2.0)/len(diffi)))

ax4.scatter(ggg, rrr, color='grey', alpha=0.6)

r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])

text = (#'f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
           '\nCorr: ' + str(round(r_value,3)) + r'$\pm$ '+  str(round(std_err,3))+
        '\nbias: ' +  str(round(bias,3))#+
        #'\nrmse: ' +  str(round(rmse,3))+
        #'\nCorrS: ' +  str(round(r_value_s,3))
        )

ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=15,
                horizontalalignment='left', verticalalignment='top')
from scipy import stats, linspace

t1 = linspace(0,50,50)
plt.plot(t1,t1,'k--')
plt.plot(t1,t1*slope + intercept,label='BoXPol', color='green', lw=1.5)
plt.plot(t1*slope2 + intercept2,t1,label='GPM DPR', color='blue', lw=1.5)

plt.legend(loc='lower right', fontsize=15, scatterpoints= 1, numpoints=1,
               shadow=True, title='lin. Regression. (Ref.)')

plt.xlim(0,50)
plt.ylim(0,50)


plt.xlabel('GPM DPR Reflectivity (dBZ)',fontsize=ff)
plt.ylabel('BoXPol Reflectivity (dBZ)',fontsize=ff)
plt.xticks(fontsize=ff)
plt.yticks(fontsize=ff)
plt.grid(color='r')
#plt.title('Scatterplot BoXPol and DPR \n 2014-10-07--02:37:44',fontsize=ff) #RW Product Polar Stereo




plt.tight_layout()
#plt.savefig('/home/velibor/shkgpm/plot/gpm_dpr_radolan_v2_'+ZP + '.png' )
#plt.close()
plt.show()



A = R

r = np.arange(0, A.shape[1])
az = np.arange(0, A.shape[0])
# mask data array for better presentation
mask_ind = np.where(A <= np.nanmin(A))
#data[mask_ind] = np.nan
ma = np.ma.array(A, mask=np.isnan(A))

from pcc import get_miub_cmap as my_cmap

cgax, caax, paax, pm = wradlib.vis.plot_cg_ppi(ma, r, az, autoext=True,

                                           refrac=False, cmap=my_cmap())
plt.tick_params(
    axis='both',
    which='both',
    bottom='off',
    top='off',
    labelbottom='off',
    right='off',
    left='off',
    labelleft='off')
t = plt.title('BoXPol PPI')
t.set_y(1.05)
cbar = plt.gcf().colorbar(pm, pad=0.075, orientation='horizontal',shrink=0.7)
plt.text(1.0, 1.05, 'azimuth', transform=caax.transAxes, va='bottom',
    ha='right')
cbar.set_label('Reflectivity (dBZ)')
plt.tight_layout()

plt.show()