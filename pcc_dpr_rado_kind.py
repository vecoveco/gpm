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

GGG = []
RRR = []

# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 12#18#7

'''
zz = np.array([20140609, 20140610, 20140629, 20140826, 20140921, 20141007,
               20141016, 20150128, 20150227, 20150402, 20150427, 20160405,
               20160607, 20160805, 20160904, 20160917, 20161001, 20161024,
               20170113, 20170203,20170223])
'''
zz = np.array(['20140709'])
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

    pfad_gpm_g = glob.glob("/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118."+str(year)+str(m)+str(d)+"*.HDF5")[0]

    print pfad_gpm_g

    gpmdpr = h5py.File(pfad_gpm_g, 'r')
    gprof_lat = np.array(gpmdpr['NS']['Latitude'])
    gprof_lon = np.array(gpmdpr['NS']['Longitude'])

    gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])



    gprof_pia = np.array(gpmdpr['NS']['SLV']['piaFinal'])

    gprof_pp[gprof_pp==-9999.9]= np.nan
    gprof_pia[gprof_pia==-9999.9]= np.nan

    print gprof_pp.shape, gprof_pia.shape
    #gprof_pp = gprof_pp + wradlib.trafo.idecibel(gprof_pia)
    #gprof_pp = gprof_pp + gprof_pia

    gprof_typ=np.array(gpmdpr['NS']['CSF']['typePrecip'])
    gprof_typ[gprof_typ==-9999] = 0
    gprof_typ[gprof_typ==-1111] = 0
    gt = gprof_typ/100000000
    gt[gt!=1]=0

    gprof_pp = gprof_pp * gt



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
    from satlib import read_rado
    #x1,y1,r1 = read_rado('201502270820')
    #rwdata = (rwdata+r1)/2



    ## Cut the GPM Swath
    ## ------------------


    blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp, eu=0)

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

    result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=400)

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

    ## Dynamischer Threshold
    THref = np.nanmax([np.nanmin(rrr),np.nanmin(ggg)])

    ## Nur Niederschlagsrelevante
    rrr[rrr < THref]=np.nan
    ggg[ggg < THref]=np.nan


    ################################################################Swap!
    #rrr, ggg = ggg, rrr

    ff = 15
    cc = 0.5
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------

    pm1 = plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)

    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("Reflectivity [dBZ]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plot_borders(ax1)
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
                         cmap=my_cmap, vmin=0.01, vmax=50, zorder=2)
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

    pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                         cmap=my_cmap, vmin=0.01, vmax=50,zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("Reflectivity [dBZ]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plt.title('RADOLAN Reflectivity Interpoliert: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
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
    slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
    line = slope * ggg +intercept

    from pcc import skill_score
    SS = skill_score(ggg,rrr,th=TH_ref)

    ax4.scatter(ggg, rrr, label='Reflectivity [dBZ]', color='grey', alpha=0.6)

    r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])

    text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
               '\nCorr: ' + str(round(r_value,3)) + r'$\pm$: '+  str(round(std_err,3))+
            '\nN: '+ str(int(SS['N']))+
            '\nHit: ' + str(SS['H'])+
            '\nMiss: ' + str(SS['M'])+
            '\nFalse: ' + str(SS['F'])+
            '\nCnegative: ' + str(SS['C'])+
            '\nHR: ' + str(round(SS['HR'],3))+
            '\nPOD: ' + str(round(SS['POD'],3))+
            '\nFAR: ' + str(round(SS['FAR'],3))+
            '\nBID: ' + str(round(SS['BID'],3))+
            '\nHSS: ' + str(round(SS['HSS'],3))+
            '\nBias: '+ str(round(SS['bias'],3))+
            '\nRMSE: '+ str(round(SS['RMSE'],3))+
            '\nCorrS:' +  str(round(r_value_s,3))
            )

    ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=10,
                    horizontalalignment='left', verticalalignment='top')

    t1 = linspace(0,50,50)
    plt.plot(t1,t1,'k-')
    plt.plot(t1, t1*slope + intercept, 'r-', lw=3 ,label='Regression')
    plt.plot(t1, t1*slope + (intercept+5), 'r-.', lw=1.5 ,label=r'Reg $\pm$ 5 mdBZ')
    plt.plot(t1, t1*slope + (intercept-5), 'r-.', lw=1.5 )
    plt.plot(np.nanmean(ggg),np.nanmean(rrr), 'ob', lw = 4,label='Mean')
    plt.plot(np.nanmedian(ggg),np.nanmedian(rrr), 'vb', lw = 4,label='Median')

    import matplotlib as mpl
    mean = [ np.nanmean(ggg),np.nanmean(rrr)]
    width = np.nanstd(ggg)
    height = np.nanstd(rrr)
    angle = 0
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                              angle=180+angle, color='blue', alpha=0.8,
                              fill=False, ls='--', label='Std')
    ax4.add_patch(ell)

    plt.legend(loc='lower right', fontsize=10, scatterpoints= 1, numpoints=1, shadow=True)

    plt.xlim(0,50)
    plt.ylim(0,50)


    plt.xlabel('GPM DPR Reflectivity [dBZ]',fontsize=ff)
    plt.ylabel('RADOLAN Reflectivity [dBZ]',fontsize=ff)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')



    plt.tight_layout()
    #plt.savefig('/home/velibor/shkgpm/plot/gpm_dpr_radolan_v2_'+ZP + '.png' )
    #plt.close()
    plt.show()

    from satlib import cp_dist
    cp_dist(ggg[maske],rrr[maske])








