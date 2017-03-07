"""

Einlesen und darstellen von GPM IMERG und Radolan Dateien


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


ZP = '20141007'#'20140609'
#year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

ye = ZP[2:4]

R_all = []
G_all_HQ = []
G_all_IR = []
G_all_MI = []
## Read GPM Data
## -------------
pfad2 = ('/home/velibor/shkgpm/data/imerg/'+str(year)+str(m)+str(d)+'/*.HDF5')
pfad_gpm = sorted(glob.glob(pfad2))

for jjj in range(len(pfad_gpm)):
    pfad_gpm_g = pfad_gpm[jjj]
    imerg_zeit = pfad_gpm_g[62:86]
    print imerg_zeit
    gpmimerg = h5py.File(pfad_gpm_g, 'r')
    gpm_lat=np.array(gpmimerg['Grid']['lat'])
    gpm_lon=np.array(gpmimerg['Grid']['lon'])

    gpm_time =np.array(gpmimerg['Grid']['HQobservationTime'], dtype=float)
    gpm_time[gpm_time==-9999]= 0

    #gpm_pp=np.array(gpmimerg['Grid']['HQprecipitation'])
    gpm_pp=np.array(gpmimerg['Grid']['precipitationCal'])
    gpm_pp_ir=np.array(gpmimerg['Grid']['IRprecipitation'])
    gpm_pp_mi=np.array(gpmimerg['Grid']['HQprecipitation'])

    gpm_pp[gpm_pp==-9999.9]= np.nan
    gpm_pp_ir[gpm_pp_ir==-9999.9]= np.nan
    gpm_pp_mi[gpm_pp_mi==-9999.9]= np.nan
    gpm_pp = np.transpose(gpm_pp)
    gpm_pp_ir = np.transpose(gpm_pp_ir)
    gpm_pp_mi = np.transpose(gpm_pp_mi)

    # Cut -------------------------------------------------------------------------

    #lim1 [xmin, xmax, ymin,ymax]
    #limit= np.array([-140, -90, -10, 30])#pazifik
    #limit= np.array([65, 100, 5, 40]) #indien
    #limit= np.array([5, 9, 49, 53]) #bonn
    #limit= np.array([-140, -80, -5, 30])#pazifik
    #limit= np.array([2,18,44,58]) #deutschland
    limit= np.array([2.0735617005681379,
                     15.704155593113517,
                     46.952580411190304,
                     54.896591448461479]) #RADOLAN
    #limit= np.array([-180,180,-90, 90]) #Welt
    lon_ii = np.array(np.where((limit[0] < gpm_lon) & (limit[1] > gpm_lon)))
    lat_ii = np.array(np.where((limit[2] < gpm_lat) & (limit[3] > gpm_lat)))

    gpm_lon = gpm_lon[lon_ii[0,0]:lon_ii[0,-1]]

    gpm_lat = gpm_lat[lat_ii[0,0]:lat_ii[0,-1]]

    gpm_pp = gpm_pp[lat_ii[0,0]:lat_ii[0,-1], lon_ii[0,0]:lon_ii[0,-1]]
    gpm_pp_ir = gpm_pp_ir[lat_ii[0,0]:lat_ii[0,-1], lon_ii[0,0]:lon_ii[0,-1]]
    gpm_pp_mi = gpm_pp_mi[lat_ii[0,0]:lat_ii[0,-1], lon_ii[0,0]:lon_ii[0,-1]]

    print gpm_pp.shape, gpm_pp_ir.shape, gpm_pp_mi.shape

    gpm_time = gpm_time[lat_ii[0,0]:lat_ii[0,-1], lon_ii[0,0]:lon_ii[0,-1]]

    print gpm_time

    if np.sum(gpm_time)==0:
        ht, mt = imerg_zeit[10:12], str(int(imerg_zeit[12:14])+15)
        print mt
    if np.sum(gpm_time)!=0:
        ht, mt = imerg_zeit[10:12], str(((int(np.nanmax(gpm_time))/5)*5)+15)
        print mt

    ## Read RADOLAN Data
    ## -----------------

    #TODO: Mean von 30 min RADOLAN nehmen! [lit. Oliveira, Charactersistics and Durial....]

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

    print radolan_zeit

    #Binaere Grid
    rn = rwdata.copy()
    rn[rn != -9999] = 1
    rn[rn == -9999] = 0

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]
    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5 # Einheit 1/100mm pro 5min
    #rwdata[rwdata < 0] = np.nan


    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    xx, yy = wradlib.georef.reproject(x, y, projection_target=proj_wgs , projection_source=proj_stereo)

    ## INTERLOLATION
    ## --------------
    gx, gy = np.meshgrid(gpm_lon, gpm_lat)
    gk3 = wradlib.georef.epsg_to_osr(31467)

    grid_gpm_xy = np.vstack((gx.ravel(), gy.ravel())).transpose()

    xy = np.vstack((xx.ravel(), yy.ravel())).transpose()

    mask = ~np.isnan(rwdata)
    #rwdata[rwdata < 0] = np.nan

    result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata[mask].reshape(900*900,1), wrl.ipol.Idw, nnearest=4)

    result = np.ma.masked_invalid(result)

    rrr = result.reshape(gx.shape)

    Z = wradlib.trafo.idecibel(rrr)
    rrr = wradlib.zr.z2r(Z, a=200., b=1.6)

    Z2 = wradlib.trafo.idecibel(rwdata)
    rwdata = wradlib.zr.z2r(Z2, a=200., b=1.6)

     ## Interpolation of the binary Grid
    ## ------------------------------
    res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1), wrl.ipol.Idw, nnearest=4)
    res_bin = res_bin.reshape(gx.shape)

    res_bin[res_bin!=0]= 1 # Randkorrektur

    rand_y_unten = 46.952580411190304
    rand_y_oben = 54.896591448461479
    rand_x_rechts = 15.704155593113517


    rrr[np.where(gy < rand_y_unten)] = np.nan
    rrr[np.where(gy > rand_y_oben)] = np.nan
    rrr[np.where(gx > rand_x_rechts)] = np.nan

    res_bin[np.where(gy < rand_y_unten)] = np.nan
    res_bin[np.where(gy > rand_y_oben)] = np.nan
    res_bin[np.where(gx > rand_x_rechts)] = np.nan
    res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

    gpm_pp = gpm_pp * res_bin
    gpm_pp_ir = gpm_pp_ir * res_bin
    gpm_pp_mi = gpm_pp_mi * res_bin

    #Todo: Radolan ueber 30 min mitteln

    # PLot ------------------------------------------------------------------------
    from pcc import plot_world

    fig = plt.figure(figsize=(16,12))
    ax11 = fig.add_subplot(242, aspect='equal')
    pm2 = plt.pcolormesh(gpm_lon,gpm_lat,np.ma.masked_invalid(gpm_pp),
                         vmin=0.1, vmax=10,cmap=my_cmap)
    cbar = plt.colorbar(pm2, shrink=0.75, orientation='horizontal')
    cbar.set_label("RainRate [mm/h]")
    plot_world(ax11,limit[0],limit[1],limit[2],limit[3])
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title('IMERG Multisat Gaug Calibration :\n ' + imerg_zeit)
    plt.xlim(limit[0],limit[1])
    plt.ylim(limit[2],limit[3])
    plt.grid()


    ax12 = fig.add_subplot(241, aspect='equal')
    pm2 = plt.pcolormesh(gpm_lon, gpm_lat,np.ma.masked_invalid(rrr),
                         vmin=0.1, vmax=10,cmap=my_cmap)
    cbar = plt.colorbar(pm2, shrink=0.75, orientation='horizontal')
    cbar.set_label("RainRate [mm/h]")
    plot_world(ax12,limit[0],limit[1],limit[2],limit[3])
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title('RADOLAN IMERG GRID: \n' + ZP +'--' + ht+':'+mt)
    plt.xlim(limit[0],limit[1])
    plt.ylim(limit[2],limit[3])
    plt.grid()

    ax13 = fig.add_subplot(243, aspect='equal')
    pm2 = plt.pcolormesh(gpm_lon,gpm_lat,np.ma.masked_invalid(gpm_pp_ir),
                         vmin=0.1, vmax=10,cmap=my_cmap)
    cbar = plt.colorbar(pm2, shrink=0.75, orientation='horizontal')
    cbar.set_label("RainRate [mm/h]")
    plot_world(ax13,limit[0],limit[1],limit[2],limit[3])
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title('IMERG IR : \n' + imerg_zeit)
    plt.xlim(limit[0],limit[1])
    plt.ylim(limit[2],limit[3])
    plt.grid()

    ax14 = fig.add_subplot(244, aspect='equal')
    pm2 = plt.pcolormesh(gpm_lon,gpm_lat,np.ma.masked_invalid(gpm_pp_mi),
                         vmin=0.1, vmax=10,cmap=my_cmap)
    cbar = plt.colorbar(pm2, shrink=0.75, orientation='horizontal')
    cbar.set_label("RainRate [mm/h]")
    plot_world(ax14,limit[0],limit[1],limit[2],limit[3])
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title('IMERG Microwave : \n' + imerg_zeit)
    plt.xlim(limit[0],limit[1])
    plt.ylim(limit[2],limit[3])
    plt.grid()

    ax15 = fig.add_subplot(245, aspect='equal')
    pm25 = plt.pcolormesh(xx, yy,rwdata,
                         vmin=0.1, vmax=10,cmap=my_cmap)
    cbar = plt.colorbar(pm25, shrink=0.75, orientation='horizontal')
    cbar.set_label("RainRate [mm/h]")
    plot_world(ax15,limit[0],limit[1],limit[2],limit[3])
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title('RADOLAN: \n' + ZP +'--' + ht+':'+mt)
    plt.xlim(limit[0],limit[1])
    plt.ylim(limit[2],limit[3])
    plt.grid()

    from pcc import plot_scatter
    ax16 = fig.add_subplot(246, aspect='equal')
    try:
        plot_scatter(gpm_pp, rrr)
    except:
        pass
    import matplotlib as mpl
    mean = [ np.nanmean(gpm_pp),np.nanmean(rrr)]
    width = np.nanstd(gpm_pp)
    height = np.nanstd(rrr)
    angle = 0
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                              angle=180+angle, color='blue', alpha=0.8,
                              fill=False, ls='--', label='Std')
    ax11.add_patch(ell)
    plt.xlabel(('GPM IMERG (mm/h)'))
    plt.ylabel(('RADOLAN (mm/h)'))
    plt.grid()

    ax17 = fig.add_subplot(247, aspect='equal')
    try:
        plot_scatter(gpm_pp_ir, rrr)
    except:
        pass
    import matplotlib as mpl
    mean = [ np.nanmean(gpm_pp_ir),np.nanmean(rrr)]
    width = np.nanstd(gpm_pp_ir)
    height = np.nanstd(rrr)
    angle = 0
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                              angle=180+angle, color='blue', alpha=0.8,
                              fill=False, ls='--', label='Std')
    ax11.add_patch(ell)
    plt.xlabel(('GPM IMERG IR (mm/h)'))
    plt.ylabel(('RADOLAN (mm/h)'))
    plt.grid()

    ax18 = fig.add_subplot(248, aspect='equal')
    try:
        plot_scatter(gpm_pp_mi, rrr)
    except:
        pass
    import matplotlib as mpl
    mean = [ np.nanmean(gpm_pp_mi),np.nanmean(rrr)]
    width = np.nanstd(gpm_pp_mi)
    height = np.nanstd(rrr)
    angle = 0
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                              angle=180+angle, color='blue', alpha=0.8,
                              fill=False, ls='--', label='Std')
    ax11.add_patch(ell)
    plt.xlabel(('GPM IMERG MI (mm/h)'))
    plt.ylabel(('RADOLAN (mm/h)'))
    plt.grid()

    plt.tight_layout()
    plt.savefig('/home/velibor/shkgpm/plot/imerg/test_imerg_'+imerg_zeit + '.png' )
    #plt.show()
    plt.close()

    R_all.append(rrr.reshape(rrr.shape[0]*rrr.shape[1]))
    G_all_HQ.append(gpm_pp.reshape(gpm_pp.shape[0]*gpm_pp.shape[1]))
    G_all_IR.append(gpm_pp_ir.reshape(gpm_pp_ir.shape[0]*gpm_pp_ir.shape[1]))
    #G_all_MI.append(gpm_pp_mi.reshape(gpm_pp_mi.shape[0]*gpm_pp_mi.shape[1]))


#print R_all.shape, G_all_HQ.shape, G_all_IR.shape, G_all_MI.shape
G_all_HQ = np.concatenate(G_all_HQ,axis=0)
G_all_IR = np.concatenate(G_all_IR,axis=0)
#G_all_MI = np.concatenate(G_all_MI,axis=0)
R_all = np.concatenate(R_all,axis=0)

fig = plt.figure(figsize=(12,12))
ax111 = fig.add_subplot(121, aspect='equal')
plot_scatter(G_all_HQ, R_all)
import matplotlib as mpl
mean = [ np.nanmean(G_all_HQ),np.nanmean(R_all)]
width = np.nanstd(G_all_HQ)
height = np.nanstd(R_all)
angle = 0
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                          angle=180+angle, color='blue', alpha=0.8,
                          fill=False, ls='--', label='Std')
ax111.add_patch(ell)
plt.xlabel(('GPM IMERG HQ (mm/h)'))
plt.ylabel(('RADOLAN (mm/h)'))
plt.grid()

ax121 = fig.add_subplot(122, aspect='equal')
plot_scatter(G_all_IR, R_all)
import matplotlib as mpl
mean = [ np.nanmean(G_all_IR),np.nanmean(R_all)]
width = np.nanstd(G_all_IR)
height = np.nanstd(R_all)
angle = 0
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                          angle=180+angle, color='blue', alpha=0.8,
                          fill=False, ls='--', label='Std')
ax121.add_patch(ell)
plt.xlabel(('GPM IMERG IR (mm/h)'))
plt.ylabel(('RADOLAN (mm/h)'))
plt.grid()



plt.savefig('/home/velibor/shkgpm/plot/imerg/all_gpm_imerg_radolan_'+ZP + '.png' )
#plt.show()
plt.close()
