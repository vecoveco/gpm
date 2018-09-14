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
from satlib import read_rado
from wradlib.trafo import idecibel
from wradlib.trafo import decibel


from satlib import writeskill2hdf as w2h


# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 0.1
scc = ['NS', 'HS', 'MS']

pfad = ('/automount/ags/velibor/gpmdata/dprV7/*.HDF5')
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

        gpmdpr_z = np.array(gpmdpr[sc]['SLV']['zFactorCorrectedNearSurface'])
        gpmdpr_z[gpmdpr_z <= 0] = np.nan
        gpmdpr_pp = np.array(gpmdpr[sc]['SLV']['precipRateNearSurface'])
        gpmdpr_pp[gpmdpr_pp <= 0] = np.nan
        #######################################################################
        dpr_pp=np.array(gpmdpr[sc]['SLV']['zFactorCorrected'])
        dpr_pp[dpr_pp < 0] = np.nan
        dpr_bbh=np.array(gpmdpr[sc]['CSF']['heightBB'], dtype=float)
        dpr_bbh[dpr_bbh < 0] = np.nan
        dpr_bbw=np.array(gpmdpr[sc]['CSF']['widthBB'], dtype=float)
        dpr_bbw[dpr_bbw < 0] = np.nan
        #phaseNearSurface
        dpr_phase=np.array(gpmdpr[sc]['SLV']['phaseNearSurface'], dtype=int)
        dpr_phase = dpr_phase/100
        #dpr_phase[dpr_phase==255]=np.nan
        #typePrecip
        dpr_typ = np.array(gpmdpr[sc]['CSF']['typePrecip'], dtype=float)
        #StromTop
        dpr_top = np.array(gpmdpr[sc]['PRE']['heightStormTop'], dtype=float)
        dpr_top[dpr_top < 0] = np.nan
        ########################################################################
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

            elif mt == '60':
                mt = '55'

            ################################################# Read RADOLAN Data
            r_pro = 'ry'

            pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
                    str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
                    str(ye)+str(m)+ str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

            pfad_radolan = pfad[:-3]

            r_pro2 = 'rx'

            pfad2 = ('/automount/radar/dwd/'+ r_pro2 +'/'+str(year)+'/'+str(year)+'-'+
                    str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro2+'_10000-'+
                    str(ye)+str(m)+ str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

            pfad_radolan2 = pfad2[:-3]


            try:

                rw_filename = wradlib.util.get_wradlib_data_file(glob.glob(pfad_radolan+'*')[0])
                rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

                rw_filename2 = wradlib.util.get_wradlib_data_file(glob.glob(pfad_radolan2+'*')[0])
                rwdata2, rwattrs2 = wradlib.io.read_RADOLAN_composite(rw_filename2)

                radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")
                #Binaere Grid
                rn = rwdata.copy()
                rn[rn != -9999] = 1
                rn[rn == -9999] = 0

                radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
                x = radolan_grid_xy[:,:,0]
                y = radolan_grid_xy[:,:,1]
                #################################### RADOLAN RX und RY Einlesen
                rwdata2 = np.ma.masked_equal(rwdata2, -9999) / 2 - 32.5
                rwdata = np.ma.masked_equal(rwdata, -9999) * 8
                #rwdata[rwdata < 0] = np.nan


                ############################################# Cut the GPM Swath

                blon, blat, gpmdpr_z_b = cut_the_swath(gprof_lon,gprof_lat,gpmdpr_z,eu=0)
                blon, blat, gpmdpr_pp_b = cut_the_swath(gprof_lon,gprof_lat,gpmdpr_pp,eu=0)
                blon, blat, dpr_pp_b = cut_the_swath(gprof_lon,gprof_lat,dpr_pp,eu=0)
                blon, blat, dpr_bbh_b = cut_the_swath(gprof_lon,gprof_lat,dpr_bbh,eu=0)
                blon, blat, dpr_bbw_b = cut_the_swath(gprof_lon,gprof_lat,dpr_bbw,eu=0)
                blon, blat, dpr_phase_b = cut_the_swath(gprof_lon,gprof_lat,dpr_phase,eu=0)

                blon, blat, dpr_typ_b = cut_the_swath(gprof_lon,gprof_lat,dpr_typ,eu=0)
                blon, blat, dpr_top_b = cut_the_swath(gprof_lon,gprof_lat,dpr_top,eu=0)

                proj_stereo = wrl.georef.create_osr("dwd-radolan")
                proj_wgs = osr.SpatialReference()
                proj_wgs.ImportFromEPSG(4326)

                gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
                grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


                ####################################### Neu
                #rwdata2[rwdata2 <= 15] = -9999


                rwdata2 = idecibel(rwdata2)



                ############################################## INTERLOLATION RY

                gk3 = wradlib.georef.epsg_to_osr(31467)

                grid_gpm_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()

                xy = np.vstack((x.ravel(), y.ravel())).transpose()

                #mask = ~np.isnan(rwdata)

                result = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata.reshape(900*900,1), wrl.ipol.Idw, nnearest=8)

                result = np.ma.masked_invalid(result)

                rrr = result.reshape(gpm_x.shape)

                ############################################## INTERLOLATION RX

                #mask2 = ~np.isnan(rwdata2)

                result2 = wrl.ipol.interpolate(xy, grid_gpm_xy, rwdata2.reshape(900*900,1), wrl.ipol.Idw, nnearest=8)

                result2 = np.ma.masked_invalid(result2)

                rrr2 = result2.reshape(gpm_x.shape)

                rrr2 = decibel(rrr2)
                rwdata2 = decibel(rwdata2)


                ## Interpolation of the binary Grid

                res_bin = wrl.ipol.interpolate(xy, grid_gpm_xy, rn.reshape(900*900,1), wrl.ipol.Idw, nnearest=25)
                res_bin = res_bin.reshape(gpm_x.shape)

                res_bin[res_bin != 0] = 1  # Randkorrektur

                rand_y_unten = -4658.6447242655722
                rand_y_oben = -3759.6447242655722
                rand_x_rechts = 375.5378330781441


                rrr[np.where(gpm_y < rand_y_unten)] = np.nan
                rrr[np.where(gpm_y > rand_y_oben)] = np.nan
                rrr[np.where(gpm_x > rand_x_rechts)] = np.nan

                rrr2[np.where(gpm_y < rand_y_unten)] = np.nan
                rrr2[np.where(gpm_y > rand_y_oben)] = np.nan
                rrr2[np.where(gpm_x > rand_x_rechts)] = np.nan

                res_bin[np.where(gpm_y < rand_y_unten)] = np.nan
                res_bin[np.where(gpm_y > rand_y_oben)] = np.nan
                res_bin[np.where(gpm_x > rand_x_rechts)] = np.nan
                res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

                ggg_pp = np.zeros((dpr_pp_b.shape[0],dpr_pp_b.shape[1],dpr_pp_b.shape[2]))
                for jx in range(dpr_pp_b.shape[2]):
                    ggg_pp[:,:,jx] = dpr_pp_b[:,:,jx] * res_bin

                ggg_bbh = dpr_bbh_b * res_bin
                ggg_bbw = dpr_bbw_b * res_bin
                ggg_phase = dpr_phase_b * res_bin
                ggg_typ = dpr_typ_b * res_bin
                ggg_pp = gpmdpr_pp_b * res_bin
                ggg_z = gpmdpr_z_b * res_bin
                ggg_top = dpr_top_b * res_bin

                #name = '/automount/ags/velibor/gpmdata/dumpdata/gpm_dpr_rado_all/' \
                #       'dprrado_'+sc+'/dprrado_'+ sc + str(gpm_zeit)
                name = '/automount/ags/velibor/gpmdata/dumpdata/gpm_dpr_rado_all_int25/' \
                       'dprrado_'+sc+'/dprrado_'+ sc + str(gpm_zeit)

                w2h(name,gpm_x, gpm_y, ggg_pp, ggg_bbh, ggg_bbw, ggg_typ,
                    ggg_phase,ggg_pp,ggg_z,ggg_top, rrr, rrr2)

                print '/automount/ags/velibor/gpmdata/dumpdataV7/SKILL_'+sc+ str(gpm_zeit)


                del(gprof_lat,
                        gprof_lon, gpmdpr_pp, res_bin, rrr, ggg_pp, ggg_bbh, ggg_bbw,
                    ggg_phase, ggg_typ, rwdata, x, y,
                        gpm_x, gpm_y, gpm_time, xy,grid_gpm_xy, grid_xy,
                        rrr2,  rn, rwattrs, result, pfad,
                        pfad_radolan, ht, m, d, ye ,mt, year, dpr_pp, dpr_bbh,
                    dpr_bbw, dpr_phase, dpr_typ,dpr_pp_b,dpr_bbh_b,dpr_bbw_b,dpr_phase_b,
                    dpr_typ_b, ggg_z,gpmdpr_pp_b,gpmdpr_z_b, r_pro2, pfad2, pfad_radolan2,
                    rwattrs2, rwdata2, rw_filename, rw_filename2, result2,name, dpr_top,
                    dpr_top_b, ggg_top)

            except:
                print 'Datei '+pfad_radolan+' nicht gefunden!'


from pcc import melde_dich
melde_dich('Das Program ist fertig!')



