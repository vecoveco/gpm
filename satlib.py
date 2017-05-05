"""

Lib zum bearbeiten von Satelliten Daten

"""


import h5py
import numpy as np
import glob
import wradlib
import datetime as dt
from osgeo import osr


def read_gprof(gprof_pfad):
    """
    Function:
        Reading gprof Koordinates and Data

    Input:
        gprof_pfad ::: Pfad zur Datei

    Output:
        Latitude, Longitude, GPROF Rainrate (mm/h)
    """
    gprof = h5py.File(gprof_pfad, 'r')
    gpmgmi_S1=gprof['S1']
    gprof_lat=np.array(gpmgmi_S1['Latitude'])
    gprof_lon=np.array(gpmgmi_S1['Longitude'])
    gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
    gprof_pp[gprof_pp<=0] = np.nan

    return gprof_lat, gprof_lon, gprof_pp


def read_dpr(dpr_pfad,scan):
    """
    Function:
        Reading DPR Data

    Input:
        dpr_pfad ::: pfad zur Datei

        scan ::: 'NS','HS','MS'

    Output:
        Latitude, Longitude, DPR Precipitation
    """
    dpr = h5py.File(dpr_pfad, 'r')
    dpr_lat=np.array(dpr[scan]['Latitude'])
    dpr_lon=np.array(dpr[scan]['Longitude'])
    dpr_pp=np.array(dpr[scan]['SLV']['zFactorCorrectedNearSurface'])
    dpr_time = dpr['NS']['ScanTime']
    #dpr_pp[dpr_pp<=0] = np.nan

    return dpr_lat, dpr_lon, dpr_pp, dpr_time


def read_rado(zeitstempel, r_pro):
    """
    Function:
        Reading RADOLAN Coordinates, Data and BinGrid

    Input:
        zeitstempel :::'YYYYMMDDhhmmss'

        r_pro ::: rx (dbz), ry (rr/h), rz(rr/h)

    Output:
        X, Y, Reflectivity/Rainrate, BinGrid
    """

    ZP = zeitstempel

    # Zeitstempel nach YYYYMMDDhhmmss
    year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
    ht, mt, st = ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]

    mt = str(int(round(float(mt)/5.0)*5.0))
    print ht, mt
    if mt == '0':
        mt = '00'
    elif mt == '5':
        mt = '05'
    elif mt =='60':
        mt = '55'



    pfad = ('/automount/radar/dwd/'+r_pro+'/'+str(year)+'/'+str(year)+'-'+str(m)+
            '/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+
            str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

    pfad_radolan = pfad[:-3]

    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    rn = rwdata.copy()
    rn[rn != -9999] = 1
    rn[rn == -9999] = 0

    if r_pro == 'rx':

        rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5
    else:
        rwdata = np.ma.masked_equal(rwdata, -9999)*8


    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]

    return x, y, rwdata, rn



def cut_the_swath(gprof_lon, gprof_lat, gprof_pp,eu):
    """
    Function:
        Cutting the ScanSwath of GPM

    Input:
        gprof_lat ::: GPM Latitude
        gprof_lon ::: GPM Longitude
        gprof_pp ::: GPM Product
        eu ::: Region Radolan(0),
               Radolan EU(1) and
               BoxPol(2)

    Output:
        blon ::: New Longitude
        blat ::: New Latitude
        gprof_pp_b ::: Cutted GPM Product

    """

    # Rand bestimmt nach Radolan Eckpunkten
    if eu==0:
        bonn_lat1 = 46.952580411190304
        bonn_lat2 = 54.896591448461479
        bonn_lon1 = 2.0735617005681379
        bonn_lon2 = 15.704155593113517
    # Rand bestimmt nach Radolan EU Eckpunkten
    if eu==1:
        bonn_lat1 = 43.874791353919626
        bonn_lat2 = 57.100558552767012
        bonn_lon1 = -0.86239071542899981
        bonn_lon2 = 21.680045338521435
    # Rand bestimmt nach BoxPol Eckpunkten
    if eu==2:
        bonn_lat1 = 49.9400
        bonn_lat2 = 51.3500
        bonn_lon1 = 6.40000
        bonn_lon2 = 8.10000

    ilat= np.where((gprof_lat>bonn_lat1) & (gprof_lat<bonn_lat2))
    ilon= np.where((gprof_lon>bonn_lon1) & (gprof_lon<bonn_lon2))

    latstart = ilat[0][0]
    latend = ilat[0][-1]


    alon = gprof_lon[latstart:latend]
    alat = gprof_lat[latstart:latend]
    gprof_pp_a = gprof_pp[latstart:latend]

    ailat= np.where((alat>bonn_lat1) & (alat<bonn_lat2))
    ailon= np.where((alon>bonn_lon1) & (alon<bonn_lon2))
    alonstart = ailon[0][0]
    alonend = ailon[0][-1]


    blon = alon[alonstart:alonend]
    blat = alat[alonstart:alonend]
    gprof_pp_b = gprof_pp_a[alonstart:alonend]

    return blon, blat, gprof_pp_b


def proj_gpm2radolan(gpm_longitude, gpm_latitude):
    """
    Function:
        Projection of cutted GPM Swath Coordinates on RadolanGrid

    Input:
        gpm_longitude ::: GPM Longitude
        gpr_latitude ::: GPM Latitude

    Output:
        gpm_x ::: GPM X Coordinate
        gpm_y ::: GPM Y Coordinate
    """
    proj_stereo = wradlib.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    gpm_x, gpm_y = wradlib.georef.reproject(gpm_longitude, gpm_latitude,
                                            projection_target=proj_stereo ,
                                            projection_source=proj_wgs)
    return gpm_x, gpm_y



def ipol_rad2gpm(radolan_x, radolan_y, gpm_x, gpm_y, radolan_data):
    """
    Function:
        Interpolation on the GPM Grid

    Input:
        radolan_x ::: RADOLAN X Coordinate
        radolan_y ::: RADOLAN Y Coordinate
        gpm_x ::: GPM X Coordinate
        gpm_y ::: GPM Y Coordinate
        radolan_daten ::: Data RADOLAN

    Output:
        rrr ::: interpolated RADOLAN GRID

    """

    grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

    xy = np.vstack((radolan_x.ravel(), radolan_y.ravel())).transpose()

    mask = ~np.isnan(radolan_data)

    result = wradlib.ipol.interpolate(xy, grid_gpm_xy, radolan_data.reshape(
        radolan_data.shape[0]*radolan_data.shape[1],1), wradlib.ipol.Idw,
                                      nnearest=4)

    result = np.ma.masked_invalid(result)

    rrr = result.reshape(gpm_x.shape)

    return rrr




"""

!!!WORK IN PROGRESS!!!

"""

def read_rado_pm5(z, z2):

    zzz = np.nansum(np.dstack((z,z2)),2)/2

    return zzz

def read_corra():
    pass
    return 'unfertig'

def read_Tb():
    pass
    return 'unfertig'


def read_IMERG():
    pass
    return 'unfertig'


def get_time_of_gpm(gpm_lon, gpm_lat, gpm_time):
    """
    Funciton:
        Returns Time where GPM is over the Mid of RADOLAN Grid

    Input:
        gpm_lon ::: GPM Longitude
        gpm_lat ::: GPM Latitude
        gpm_time :: GPM Times

    Output:
        gpm_dt ::: GPM Time

    """
    #Todo: Verbessern
    mitte = gpm_lon.shape[1]/2 # midel swath
    ii = np.where(((gpm_lon[:,mitte]<15) & (gpm_lon[:,mitte]>2)) & ((gpm_lat[:,mitte]<54) & (gpm_lat[:,mitte] > 46)))
    gpm_year = int(np.nanmedian(np.array(gpm_time['Year'])[ii]))
    gpm_month = int(np.nanmedian(np.array(gpm_time['Month'])[ii]))
    gpm_day = int(np.nanmedian(np.array(gpm_time['DayOfMonth'])[ii]))
    gpm_hour = int(np.nanmedian(np.array(gpm_time['Hour'])[ii]))
    gpm_min = int(np.nanmedian(np.array(gpm_time['Minute'])[ii]))
    gpm_sek = int(np.nanmedian(np.array(gpm_time['Second'])[ii]))
    gpm_dt = dt.datetime(gpm_year,gpm_month, gpm_day, gpm_hour, gpm_min, gpm_sek).strftime("%Y%m%d%H%M%S")
    return gpm_dt

'''
x1,y1,z1 = read_rado("201410070235")
x2,y2,z2 = read_rado("201410070240")
z1[z1<0.2]=np.nan
z2[z2<0.2]=np.nan
from scipy import stats, linspace
maske = ~np.isnan(z1) & ~np.isnan(z2)
slope, intercept, r_value, p_value, std_err = stats.linregress(z1[maske], z2[maske])
plt.subplot(2,2,1)
plt.pcolormesh(x1,y1,z1, vmin=0.1, vmax=10)
plt.subplot(2,2,2)
plt.pcolormesh(x2,y2,z2, vmin=0.1, vmax=10)
plt.subplot(2,2,3)
plt.pcolormesh(x1,y2,z1-z2, vmin=0.1, vmax=10)
plt.subplot(2,2,4)
plt.scatter(z1,z2)
plt.title(str(r_value))
plt.show()
'''

def boxpol_pos():
    """

    Koordinaten des Bonner Radar BoXPol
    geografische Laenge 	7.071663 Ost
    geografische Breite 	50.73052 Nord
    Hohe uber NN 	99.5 m
    """

    pos_boxpol = {'lat_ppi' : 50.730519999999999,
                  'lon_ppi' : 7.071663,
                  'gky_ppi' : -4235.233235191105,
                  'gkx_ppi' : -216.64772430049572}

    return pos_boxpol


def cp_dist(data1, data2):
    """
    Function:
        Plot of PDF und CDF of two vectors

    Input:
        data1, data2 ::: Input Data

    Output:
        plot of PDF and CDF

    """
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots()
    counts1, bins1, patches1 = plt.hist(data1, alpha=0.4, color='green')
    counts2, bins2, patches2 = plt.hist(data2, alpha=0.4, color='blue')

    bin_centers1 = np.mean(zip(bins1[:-1], bins1[1:]), axis=1)
    bin_centers2 = np.mean(zip(bins2[:-1], bins2[1:]), axis=1)

    ax.plot(bin_centers1, counts1.cumsum(), 'go-')
    ax.plot(bin_centers2, counts2.cumsum(), 'bo-')


    plt.show()

