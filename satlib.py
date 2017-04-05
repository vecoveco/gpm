"""

Lib zum bearbeiten von Satelliten Daten

"""


import h5py
import numpy as np
import glob
import wradlib
import datetime as dt


def read_gprof(gprof_pfad):
    gprof = h5py.File(gprof_pfad, 'r')
    gpmgmi_S1=gprof['S1']
    gprof_lat=np.array(gpmgmi_S1['Latitude'])
    gprof_lon=np.array(gpmgmi_S1['Longitude'])
    gprof_pp=np.array(gpmgmi_S1['surfacePrecipitation'])
    gprof_pp[gprof_pp<=0] = np.nan

    return gprof_lat, gprof_lon, gprof_pp


def read_dpr(dpr_pfad):

    dpr = h5py.File(dpr_pfad, 'r')
    dpr_lat=np.array(dpr['NS']['Latitude'])
    dpr_lon=np.array(dpr['NS']['Longitude'])
    dpr_pp=np.array(dpr['NS']['SLV']['precipRateNearSurface'])
    dpr_pp[dpr_pp<=0] = np.nan

    return dpr_lat, dpr_lon, dpr_pp

def read_rado(zeitstempel):
    #ZP = '20150110220500'
    ZP = zeitstempel

    # Zeitstempel nach YYYYMMDDhhmmss
    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]

    r_pro = 'rx'

    pfad = ('/automount/radar/dwd/'+r_pro+'/'+str(year)+'/'+str(year)+'-'+str(m)+
            '/'+str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+
            str(m)+str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

    pfad_radolan = pfad[:-3]

    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    #rwdata = np.ma.masked_equal(rwdata, -9999)*8
    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]

    return x, y, rwdata


def read_rado_pm5(z, z2):

    zzz = np.nansum(np.dstack((z,z2)),2)/2

    return zzz

def read_corra():
    pass
    return 'unfertig'

def read_Tb():
    pass
    return 'unfertig'


def time_of_dpr(gpm_lon, gpm_lat, gpm_time):
    ii = np.where(((gpm_lon[:,:]<15) & (gpm_lon[:,:]>2)) & ((gpm_lat[:,:]<54) & (gpm_lat[:,:] > 46)))
    gpm_year = int(np.median(np.array(gpm_time['Year'])[ii]))
    gpm_month = int(np.median(np.array(gpm_time['Month'])[ii]))
    gpm_day = int(np.median(np.array(gpm_time['DayOfMonth'])[ii]))
    gpm_hour = int(np.median(np.array(gpm_time['Hour'])[ii]))
    gpm_min = int(np.median(np.array(gpm_time['Minute'])[ii]))
    gpm_sek = int(np.median(np.array(gpm_time['Second'])[ii]))
    gpm_dt = dt.datetime(gpm_year,gpm_month, gpm_day, gpm_hour, gpm_min, gpm_sek).strftime("%Y.%m.%d -- %H:%M:%S")
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



def cp_dist(data1, data2):
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

