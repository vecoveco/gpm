"""

Das Program dient der Veranschaulichung der 5 Minutigen RX Radolan Daten!
Es werden durch den Zeitstempel Deutschlandweite Niederschlags und
Reflektivitaeten dargestellt!


"""



import numpy as np
import wradlib




############################################### Zeitstempel nach YYYYMMDDhhmmss

from pcc import zeitschleife as zt

# YYYY MM DD hh mm mili
zeit = zt(2017,05,02,18,50,0,
          2017,05,04,23,55,0,
          steps=5)

r_data = np.zeros((len(zeit),13))


for ij in range(len(zeit)):

    print 'Zeitstempel ',ij,' von  ', len(zeit)

    ZP = zeit[ij]
    print 'Zeitpunkt: ',ZP

    year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    ye = ZP[2:4]

    ################################################### Read RADOLAN GK Koordinaten

    iii = 0
    r_pro = 'rx'
    pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
            str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+str(ye)+str(m)+
            str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

    pfad_radolan = pfad[:-3]


    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)


    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

    radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
    radolan_zeit_sav = rwattrs['datetime'].strftime("%Y%m%d-%H%M%S")

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)

    r_data[ij,0]= int(ZP)
    r_data[ij,1]= np.nanmean(rwdata[rwdata>0])
    r_data[ij,2]= np.nanstd(rwdata[rwdata>0])
    r_data[ij,3]= np.nanmedian(rwdata[rwdata>0])
    r_data[ij,4]= np.nanmax(rwdata[rwdata>0])
    r_data[ij,5]= np.nanmin(rwdata[rwdata>0])
    r_data[ij,6]= rwdata[rwdata>50].shape[0]
    r_data[ij,7]= rwdata[rwdata>40].shape[0]
    r_data[ij,8]= rwdata[rwdata>30].shape[0]
    r_data[ij,9]= rwdata[rwdata>20].shape[0]
    r_data[ij,10]= rwdata[rwdata>10].shape[0]
    r_data[ij,11]= rwdata[rwdata>1].shape[0]
    r_data[ij,12]= rwdata[rwdata>=np.nanmin(rwdata)].shape[0]




np.save('/automount/ags/velibor/plot/radolan/egon201701/r_'+radolan_zeit, r_data)