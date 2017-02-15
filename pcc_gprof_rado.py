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
#from pcc import get_miub_cmap
#my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap = get_my_cmap()


zz = np.array([20140921, 20141007,20140826,
               20141016, 20150128, 20150227, 20150402, 20150427, 20160405,
               20160607, 20160805, 20160904, 20160917, 20161001, 20161024,
               20170113, 20170203,])
#zz = np.array([20140921])

for i in range(len(zz)):
    ZP = str(zz[i])
    #year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

    ye = ZP[2:4]

    ## Read GPM Data
    ## -------------

    pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/gprof/*.HDF5')
    pfad_gprof = glob.glob(pfad2)
    pfad_gprof_g = pfad_gprof[0]

    gpmdprs = h5py.File(pfad_gprof_g, 'r')
    gprof_lat=np.array(gpmdprs['S1']['Latitude'])
    gprof_lon=np.array(gpmdprs['S1']['Longitude'])

    gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])

    gpm_time = gpmdprs['S1']['ScanTime']
    gpm_zeit = get_time_of_gpm(gprof_lon, gprof_lat, gpm_time)
    print gpm_zeit

    gprof_pp[gprof_pp==-9999.9]= np.nan

    ht, mt = gpm_zeit[14:16], str(int(round(float(gpm_zeit[17:19])/5.0)*5.0))
    if mt == '0':
        mt = '00'
    elif mt == '5':
        mt = '05'
    print mt

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

    radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")

    #Binaere Grid
    rn = rwdata.copy()
    rn[rn != -9999] = 1
    rn[rn == -9999] = 0

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]
    rwdata = np.ma.masked_equal(rwdata, -9999) *8 # Einheit 1/100mm pro 5min
    #rwdata[rwdata < 0] = np.nan


    ## Cut the GPM Swath
    ## ------------------

    from pcc import cut_the_swath
    blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp, eu=True)

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
    rrr[rrr <= 0.1] = np.nan
    ggg[ggg <= 0.1] = np.nan

    rrr = rrr*(ggg/ggg) # Rand vom GPROF entfernen !



    ff = 15
    cc = 0.5
    vmini=0
    maxi = np.nanmax([ggg,rrr])+2
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(221, aspect='equal')#------------------------------------

    pm1 = plt.pcolormesh(x, y, rwdata, cmap=my_cmap, vmin=vmini, vmax=10, zorder=2)

    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("RR [mm/h]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plot_borders(ax1)
    plot_radar(bonnlon, bonnlat, ax1, reproject=True)

    plt.title('RADOLAN: \n'+ radolan_zeit + ' UTC',fontsize=ff)
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
                         cmap=my_cmap, vmin=vmini, vmax=10, zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("RR [mm/h]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.title('GPM GPROF: \n'+'20' + gpm_zeit + ' UTC',fontsize=ff)
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

    ax2 = fig.add_subplot(223, aspect='equal')#------------------------------------

    pm3 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr),
                         cmap=my_cmap, vmin=vmini, vmax=10,zorder=2)
    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label("RR [mm/h]",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)

    plt.title('RADOLAN Interpoliert: \n'+ radolan_zeit + ' UTC',fontsize=ff) #RW Product Polar Stereo
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

    ax4 = fig.add_subplot(224, aspect='equal')#------------------------------------

    maske = ~np.isnan(ggg) & ~np.isnan(rrr)
    slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
    line = slope * ggg +intercept

    from pcc import skill_score
    SS = skill_score(ggg,rrr,th=0)

    plt.scatter(ggg, rrr, label='RR [mm/h]', color='grey', alpha=0.6)

    text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
               '\nCorr: ' + str(round(r_value,3)) + r'$\pm$: '+  str(round(std_err,3))+
            '\nN: '+ str(int(SS['N']))+
            '\nHit: ' + str(round(SS['H']/SS['N'],3)*100)+'%'+
            '\nMiss: ' + str(round(SS['M']/SS['N'],3)*100)+'%'+
            '\nFalse: ' + str(round(SS['F']/SS['N'],3)*100)+'%'+
            '\nCnegative: ' + str(round(SS['C']/SS['N'],3)*100)+'%'+
            '\nPOD: ' + str(round(SS['POD'],3))+
            '\nFAR: ' + str(round(SS['FAR'],3))+
            '\nBID: ' + str(round(SS['BID'],3))+
            '\nHSS: ' + str(round(SS['HSS'],3))+
            '\nBias: '+ str(round(SS['bias'],3))+
            '\nRMSE: '+ str(round(SS['RMSE'],3))
            )

    ax4.annotate(text, xy=(0.01, 0.99), xycoords='axes fraction', fontsize=10,
                    horizontalalignment='left', verticalalignment='top', color='blue')

    t1 = linspace(0,maxi,maxi)
    plt.plot(t1,t1,'k-')
    #plt.plot(t1,t1 + 5,'k-.')
    #plt.plot(t1,t1 - 5,'k-.')
    plt.plot(t1, t1*slope + intercept, 'r-', lw=3 ,label='Regression')
    plt.plot(t1, t1*slope + (intercept+5), 'r-.', lw=1.5 ,label='Regression + 5 mm/h')
    plt.plot(t1, t1*slope + (intercept-5), 'r-.', lw=1.5 ,label='Regression - 5 mm/h')
    plt.plot(np.nanmean(ggg),np.nanmean(rrr), 'ob', lw = 4,label='Mean')
    plt.plot(np.nanmedian(ggg),np.nanmedian(rrr), 'vb', lw = 4,label='Median')


    import matplotlib as mpl
    mean = [ np.nanmean(ggg),np.nanmean(rrr)]
    width = np.nanstd(ggg)
    height = np.nanstd(rrr)
    print width, height
    angle = 0
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                              angle=180+angle, color='blue', alpha=0.8, fill=False, ls='--')
    ax4.add_patch(ell)

    plt.xlim(0,maxi)
    plt.ylim(0,maxi)

    plt.legend(loc='lower right', fontsize=10, scatterpoints= 1, numpoints=1, shadow=True)

    plt.xlabel('GPM GPROF RR [mm/h]',fontsize=ff)
    plt.ylabel('RADOLAN RR [mm/h]',fontsize=ff)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)
    plt.grid(color='r')


    plt.tight_layout()
    plt.savefig('/home/velibor/shkgpm/plot/gpm_gprof_radolan_'+ZP + '.png' )
    plt.close()





plt.subplot(1,2,1)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='r',lw=3)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='r',lw=3)
plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr*(ggg/ggg)), cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)

plt.subplot(1,2,2)
plt.plot(gpm_x[:,0],gpm_y[:,0], color='r',lw=3)
plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='r',lw=3)
plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(rrr), cmap=my_cmap, vmin=0.01, vmax=10, zorder=2)
plt.show()
