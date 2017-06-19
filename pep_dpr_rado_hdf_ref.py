"""

Einlesen und darstellen von GPM DPR und Radolan Dateien

Radolanpfad:

"""


import h5py
import numpy as np
import wradlib
import glob
import wradlib as wrl
from osgeo import osr
from pcc import get_time_of_gpm
from pcc import cut_the_swath


from satlib import write2hdf as w2h

import time
sekunden = 60*60*5
time.sleep(sekunden)

print 'start nach ',str(sekunden), ' Sekunden.'
from pcc import melde_dich
melde_dich('Das Program REF STARTET!')

# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 0.1
scc = ['NS', 'HS', 'MS']

pfad = ('/automount/ags/velibor/gpmdata/dpr/*.HDF5')
pfad_gpm = sorted(glob.glob(pfad))

print 'Es sind ', len(pfad_gpm), ' vorhanden!'

for iii in scc:
    sc = iii

    #len(pfad_gpm)

    for i in range(0, len(pfad_gpm)):

        ## Read GPM Data

        pfad_gpm_g = pfad_gpm[i]
        print pfad_gpm_g
        gpmdpr = h5py.File(pfad_gpm_g, 'r')
        # sc = ['NS', 'HS', 'MS']

        gprof_lat = np.array(gpmdpr[sc]['Latitude'])
        gprof_lon = np.array(gpmdpr[sc]['Longitude'])

        gprof_pp = np.array(gpmdpr[sc]['SLV']['zFactorCorrectedNearSurface'])
        gprof_pp[gprof_pp<=0]= np.nan

        gpm_time = gpmdpr[sc]['ScanTime']

        try:
            gpm_zeit = get_time_of_gpm(gprof_lon, gprof_lat, gpm_time)
        except ValueError:
            pass
            print ('____________ValueError____________')
        else:
            print gpm_zeit

            ht, mt = gpm_zeit[14:16], str(int(round(float(gpm_zeit[17:19])/5.0)*5.0))
            year, ye, m, d = gpm_zeit[0:4], gpm_zeit[2:4], gpm_zeit[5:7], gpm_zeit[8:10]
            print ht, mt
            if mt == '0':
                mt = '00'
            elif mt == '5':
                mt = '05'

            elif mt =='60':
                mt = '55'

            ## Read RADOLAN Data
            r_pro = 'rx'

            pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
                    str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
                    str(ye)+str(m)+ str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

            pfad_radolan = pfad[:-3]


            try:

                rw_filename = wradlib.util.get_wradlib_data_file(glob.glob(pfad_radolan+'*')[0])

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
                #rwdata = np.ma.masked_equal(rwdata, -9999) *8
                #rwdata[rwdata < 0] = np.nan


                ## Cut the GPM Swath

                blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp,eu=0)

                proj_stereo = wrl.georef.create_osr("dwd-radolan")
                proj_wgs = osr.SpatialReference()
                proj_wgs.ImportFromEPSG(4326)

                gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
                grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


                ## INTERLOLATION

                gk3 = wradlib.georef.epsg_to_osr(31467)

                grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

                xy = np.vstack((x.ravel(), y.ravel())).transpose()

                mask = ~np.isnan(rwdata)

                result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=4)

                result = np.ma.masked_invalid(result)

                rrr = result.reshape(gpm_x.shape)


                ## Interpolation of the binary Grid

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

                w2h('/automount/ags/velibor/gpmdata/dumpdata/REF_'+sc+ str(gpm_zeit), gpm_x, gpm_y, rrr, ggg)


                del(gprof_lat,
                        gprof_lon, gprof_pp, res_bin, rrr, ggg, rwdata, x, y,
                        gpm_x, gpm_y, gpm_time, xy,grid_gpm_xy, grid_xy,
                        mask,  rn, rwattrs, result,  pfad,
                        pfad_radolan, ht, m, d, ye ,mt, year)

            except:
                print 'Datei '+pfad_radolan+' nicht gefunden!'


from pcc import melde_dich
melde_dich('Das Program REF ist fertig!')




