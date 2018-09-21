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
import scipy as sp
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
from pcc import get_my_cmap
my_cmap = get_my_cmap()

from pcc import get_my_cmap
my_cmap2 = get_my_cmap()

GGG = []
RRR = []

from pcc import get_radar_locations
from pcc import plot_radar2

def plot_all_cband(ax):
    for i in get_radar_locations().keys():


        plot_radar2(get_radar_locations()[i]['lon'],
                    get_radar_locations()[i]['lat'], ax , reproject=True, cband=True, col='black')

# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 12#18#7

'''
zz = np.array([20140609, 20140610, 20140629, 20140826, 20140921, 20141007,
               20141016, 20150128, 20150227, 20150402, 20150427, 20160405,
               20160607, 20160805, 20160904, 20160917,
               20170113])
'''
zz = np.array(['20150227'])


for i in range(len(zz)):
    ZP = str(zz[i])
    #year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

    ye = ZP[2:4]

    ## Read GPM Data
    ## -------------

    #pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
    #pfad_gpm = glob.glob(pfad2)
    #pfad_gpm_g = pfad_gpm[0]

    #pfad_gpm_g = glob.glob("/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118."+str(year)+str(m)+str(d)+"*.HDF5")[0]
    pfad_gpm_g = glob.glob("/automount/ags/velibor/gpmdata/dprV7/2A.GPM.DPR.V7-20170308."+str(year)+str(m)+str(d)+"*.HDF5")[0]

    print pfad_gpm_g

    gpmdpr = h5py.File(pfad_gpm_g, 'r')
    gprof_lat = np.array(gpmdpr['NS']['Latitude'])
    gprof_lon = np.array(gpmdpr['NS']['Longitude'])

    gprof_pp = np.array(gpmdpr['NS']['SLV']['precipRateESurface'])
    #gprof_pp = np.array(gpmdpr['NS']['SLV']['precipRateNearSurface'])
    #gprof_pp = np.array(gpmdpr['NS']['SLV']['precipRateAve24'])

    #gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])
    #gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrectedESurface'])
    #gprof_pp = np.array(gpmdpr['NS']['SLV']['precipRateNearSurface'])



    gpm_h = np.array(gpmdpr['NS']['PRE']['elevation'])
    gpm_h[gpm_h ==-9999.9]= np.nan

    gprof_pp[gprof_pp ==-9999.9]= np.nan

    #gprof_pp = gprof_pp + wradlib.trafo.idecibel(gprof_pia)
    #gprof_pp = gprof_pp + gprof_pia



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

    r_pro = 'ry'

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
    rwdata = np.ma.masked_equal(rwdata, -9999) *8#/ 2 - 32.5

    rwdata[rwdata <= 0.5] = -9999


    #print ('Min RadolaN: ',np.nanmin(rwdata))
    #rwdata = np.log10(rwdata)
    from satlib import read_rado
    #x1,y1,r1 = read_rado('201502270820')
    #rwdata = (rwdata+r1)/2
    #from wradlib.trafo import idecibel
    #from wradlib.trafo import decibel
    #rwdata = idecibel(rwdata)



    ## Cut the GPM Swath
    ## ------------------


    blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp, eu=0)
    blon1, blat1, gpm_h = cut_the_swath(gprof_lon,gprof_lat,gpm_h, eu=0)############################################################

    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
    grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


    #rwdata = np.log10(rwdata)

    ## INTERLOLATION
    ## --------------

    gk3 = wradlib.georef.epsg_to_osr(31467)

    grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

    xy = np.vstack((x.ravel(), y.ravel())).transpose()

    mask = ~np.isnan(rwdata)

    result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=8)

    result = np.ma.masked_invalid(result)

    rrr = result.reshape(gpm_x.shape)

    #rrr = decibel(rrr)
    #rwdata = decibel(rwdata)
    #rrr = 10**rrr
    #rwdata = 10**rwdata



    ## Interpolation of the binary Grid
    ## ------------------------------
    #rhmax = np.load("/automount/ags/velibor/data/radolan_dx/RY_maxHxy_testradius1500m.npy")
    #rhmin = np.load("/automount/ags/velibor/data/radolan_dx/RY_minHxy.npy")
    #ry_h = np.ma.masked_invalid(rhmax[2]-rhmin[2])

    ry_h = np.load("/automount/ags/velibor/data/radolan_dx/RY_minHxy.npy")
    ry_h = np.ma.masked_invalid(ry_h[2])
    print ('-------------------------------->', ry_h.shape)


    res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1), wrl.ipol.Idw, nnearest=8)
    res_bin = res_bin.reshape(gpm_x.shape)

    ry_hi = wrl.ipol.interpolate(xy, grid_gpm_xy, ry_h.reshape(900*900,1), wrl.ipol.Idw, nnearest=8)
    ry_hi = ry_hi.reshape(gpm_x.shape)

    res_bin[res_bin!=0]= 1 # Randkorrektur

    rand_y_unten = -4658.6447242655722
    rand_y_oben = -3759.6447242655722
    rand_x_rechts = 375.5378330781441


    ry_hi[np.where(gpm_y < rand_y_unten)] = np.nan
    ry_hi[np.where(gpm_y > rand_y_oben)] = np.nan
    ry_hi[np.where(gpm_x > rand_x_rechts)] = np.nan

    rrr[np.where(gpm_y < rand_y_unten)] = np.nan
    rrr[np.where(gpm_y > rand_y_oben)] = np.nan
    rrr[np.where(gpm_x > rand_x_rechts)] = np.nan

    res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
    res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
    res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan
    res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

    ggg = gprof_pp_b * res_bin
    ggh = gpm_h * res_bin
    ry_hi = ry_hi * res_bin

    print ('-------------------------------->', ry_hi.shape)

    ## Dynamischer Threshold
    #THref = np.nanmax([np.nanmin(rrr),np.nanmin(ggg)])

    ## Nur Niederschlagsrelevante
    #rrr[rrr < THref]=np.nan
    #ggg[ggg < THref]=np.nan

    # Normalisieren
    #rrr = rrr/np.nanmax(rrr)
    #ggg = ggg/np.nanmax(ggg)
    ggg[ggg<=0.5]=np.nan

    ################################################################Swap!
    #rrr, ggg = ggg, rrr

    ff = 15
    cc = 0.5
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------

    pm1 = plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)

    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("Reflectivity [dBZ]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plot_borders(ax1)
    plot_all_cband(ax1)
    plot_radar(bonnlon, bonnlat, ax1, reproject=True, cband=False,col='black')


    plt.title('RADOLAN Reflectivity: \n'+ radolan_zeit + ' UTC',fontsize=ff)
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

    ax2 = fig.add_subplot(222, aspect='equal')#------------------------------------

    pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ggg),
                         cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("Reflectivity [dBZ]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.title('GPM DPR Reflectivity: \n'+ gpm_zeit + ' UTC',fontsize=ff)
    plot_borders(ax2)
    plot_radar(bonnlon, bonnlat, ax2, reproject=True, cband=False,col='black')
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


    ax2 = fig.add_subplot(223, aspect='equal')#------------------------------------

    pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ry_hi),
                         cmap='jet', vmin=500, vmax=3000,zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("Hight in m",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plt.title('RADOLAN Height Interpoliert: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
    plot_borders(ax2)
    plot_radar(bonnlon, bonnlat, ax2, reproject=True, cband=False,col='black')
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

    ax4 = fig.add_subplot(224, aspect='equal')#------------------------------------


    maske = ~np.isnan(ggg) & ~np.isnan(rrr)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
    #line = slope * ggg +intercept
    r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])


    ax4.scatter(ggg, rrr,c=ry_hi,s=30,label='Reflectivity [dBZ]', color='black', alpha=0.9)
    plt.colorbar()



    t1 = linspace(0,50,50)
    plt.plot(t1,t1,'k-')

    #plt.legend(loc='lower right', fontsize=10, scatterpoints= 1, numpoints=1, shadow=True)

    plt.xlim(0,np.nanmax([rrr,ggg]))
    plt.ylim(0,np.nanmax([rrr,ggg]))

    plt.title(r_value_s)
    plt.xlabel('GPM DPR Reflectivity [dBZ]',fontsize=ff)
    plt.ylabel('RADOLAN Reflectivity [dBZ]',fontsize=ff)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')


    plt.plot(t1,t1,'k-')

    plt.tight_layout()
    #plt.savefig('/home/velibor/shkgpm/plot/gpm_dpr_radolan_v2_'+ZP + '.png' )
    #plt.close()
    plt.show()

    #from satlib import cp_dist
    #cp_dist(ggg[maske],rrr[maske])

    #from satlib import validation_plot
    #validation_plot(ggg,rrr,15)
###################################################################################################################################
    hth = 2000
    plt.scatter(ggg[abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh))< hth],
                rrr[abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh))< hth],
                c=ry_hi[abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh))< hth],
                s=30,label='Reflectivity [dBZ]', color='black', alpha=0.9)

    plt.colorbar()
    plt.xlabel('GPM DPR Reflectivity [dBZ]',fontsize=ff)
    plt.ylabel('RADOLAN Reflectivity [dBZ]',fontsize=ff)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')
    plt.xlim(0,np.nanmax([rrr,ggg]))
    plt.ylim(0,np.nanmax([rrr,ggg]))
    plt.title(np.corrcoef(ggg[abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh))< hth],
                rrr[abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh))< hth]))

    plt.show()

    ff = 15
    cc = 0.5
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------

    pm1 = plt.pcolormesh(gpm_x, gpm_y, np.ma.masked_invalid(ggh), cmap=my_cmap, zorder=2)

    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("H in m",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plot_borders(ax1)
    plot_all_cband(ax1)
    plot_radar(bonnlon, bonnlat, ax1, reproject=True, cband=False,col='black')


    plt.title('GPM Hoehen: \n'+ radolan_zeit + ' UTC',fontsize=ff)
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

    ax2 = fig.add_subplot(222, aspect='equal')#------------------------------------

    pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ry_hi),
                         cmap=my_cmap, zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("H in m",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.title('RADOLAN Hoehen: \n'+ gpm_zeit + ' UTC',fontsize=ff)
    plot_borders(ax2)
    plot_radar(bonnlon, bonnlat, ax2, reproject=True, cband=False,col='black')
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


    ax2 = fig.add_subplot(223)#------------------------------------

    pm3 = plt.scatter(np.ma.masked_invalid(ry_hi),np.ma.masked_invalid(ggh))
    plt.xlabel('RADOLAN H')
    plt.ylabel('GPM H')

    plt.xlim(0,np.nanmax(np.ma.masked_invalid(ry_hi)))
    plt.ylim(0,np.nanmax(np.ma.masked_invalid(ggh)))


    plt.grid(color='r')

    ax234 = fig.add_subplot(224, aspect='equal')
    pm234 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh),
                         cmap=my_cmap, zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("H in m",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.title('RADOLAN - GPM H Hoehen: \n'+ gpm_zeit + ' UTC',fontsize=ff)
    plot_borders(ax234)
    plot_radar(bonnlon, bonnlat, ax234, reproject=True, cband=False,col='black')
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


    plt.show()
    print np.nanmean(abs(np.ma.masked_invalid(ry_hi)-np.ma.masked_invalid(ggh)))

