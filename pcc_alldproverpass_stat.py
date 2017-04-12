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

import csv

# Ref.Threshold nach RADOLAN_Goudenhoofdt_2016
TH_ref = 12#18#7

pfad = ('/automount/ags/velibor/gpmdata/dpr/*.HDF5')
pfad_gpm = sorted(glob.glob(pfad))

print 'Es sind ', len(pfad_gpm), ' vorhanden!'

f = open('/home/velibor/shkgpm/texte/dpr_stat_rr.csv','w')
writer = csv.writer(f, dialect='excel')
writer.writerow(['time','r_value','std_err','N','H','M','F','C','HR','POD','FAR','BID','HSS','bias','RMSE','meanG','meanR','medG','medR'])


for i in range(0, len(pfad_gpm)):
    #ZP = str(zz[i])
    #year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
    #year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
    #ye = ZP[2:4]


    ## Read GPM Data
    ## -------------
    #try:
    pfad_gpm_g = pfad_gpm[i]
    print pfad_gpm_g
    gpmdpr = h5py.File(pfad_gpm_g, 'r')
    gprof_lat = np.array(gpmdpr['NS']['Latitude'])
    gprof_lon = np.array(gpmdpr['NS']['Longitude'])

    #gprof_pp = np.array(gpmdpr['NS']['SLV']['zFactorCorrectedNearSurface'])
    #gprof_pp[gprof_pp==-9999.9]= np.nan

    gprof_pp = np.array(gpmdpr['NS']['SLV']['precipRateNearSurface'])
    gprof_pp[gprof_pp==-9999.9]= np.nan




    gpm_time = gpmdpr['NS']['ScanTime']

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


        ## Cut the GPM Swath
        ## ------------------


        blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp,eu=0)

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

        ## Dynamischer Threshold
        #THref = np.nanmax([np.nanmin(rrr),np.nanmin(ggg)])
        THref = 0.1

        ## Z in R
        Z = wradlib.trafo.idecibel(rwdata)
        rwdata = wradlib.zr.z2r(Z, a=200., b=1.6)
        #rwdata = z2r2(Z)
        #z2rEnhanced

        Z2 = wradlib.trafo.idecibel(rrr)
        rrr = wradlib.zr.z2r(Z2, a=200., b=1.6)

        ## Nur Niederschlagsrelevante
        rrr[rrr < THref]=np.nan
        ggg[ggg < THref]=np.nan


        ################################################################Swap!
        #rrr, ggg = ggg, rrr


        try:
            maske = ~np.isnan(ggg) & ~np.isnan(rrr)
            slope, intercept, r_value, p_value, std_err = stats.linregress(ggg[maske], rrr[maske])
            line = slope * ggg +intercept

            from pcc import skill_score
            SS = skill_score(ggg,rrr,th=TH_ref)

            r_value_s, p_value_s = stats.spearmanr(ggg[maske],rrr[maske])

            text = ('f(x) = ' + str(round(slope,3)) + 'x + ' + str(round(intercept,3)) +
                       '\nCorr: ' + str(round(r_value,3)) + r'$\pm$ '+  str(round(std_err,3))+
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



            writer.writerow([gpm_zeit,
                            r_value, std_err,SS['N'],SS['H'],
                            SS['M'],SS['F'],SS['C'],SS['HR'],SS['POD'],
                            SS['FAR'],SS['BID'],SS['HSS'],SS['bias'],
                            SS['RMSE'],np.nanmean(ggg),np.nanmean(rrr),
                            np.nanmedian(ggg),np.nanmedian(rrr)])

            del(text, slope, intercept, r_value, p_value, std_err, line,
                width, height, ell,maske,SS)

        except:
            pass


        #plt.show()



        del(gprof_lat, gprof_lon, gprof_pp, res_bin, rrr, ggg, rwdata, x, y,
            gpm_x, gpm_y, gpm_time, xy,grid_gpm_xy, grid_xy, mask,  rn,
            rwattrs, result,  pfad, pfad_radolan, ht, m, d, ye ,mt, year )



from pcc import melde_dich
melde_dich('Das Program pcc_alldproverpass_stats_rr.py ist fertig!')


f.close()



