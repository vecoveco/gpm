#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import math
import pandas as pd
import wradlib
from scipy import stats
import matplotlib.cm as cm
my_cmap = cm.get_cmap('jet',40)
my_cmap.set_under('lightgrey')
my_cmap.set_over('darkred')
from pcc import get_miub_cmap as my_cmap
from pcc import plot_radar
from pcc import boxpol_pos
from pcc import plot_borders
import wradlib as wrl
from osgeo import osr
from satlib import ipoli_radi_stat
from satlib import corcor
Pos = boxpol_pos()
blon0, blat0 = Pos['lon_ppi'], Pos['lat_ppi']
bbx, bby = Pos['gkx_ppi'], Pos['gky_ppi']
from time import *
from pcc import cut_the_swath
from satlib import good_overpasses_dpr_boxpol as overpasses_dpr_boxpol

print ("_")



tstart = clock()

#print (overpasses_dpr_boxpol)

for ii in overpasses_dpr_boxpol.keys()[0]:
    ZP = overpasses_dpr_boxpol[ii][0],
    pfadnr = overpasses_dpr_boxpol[ii][1]
    enigma = overpasses_dpr_boxpol[ii][2]
    offset = overpasses_dpr_boxpol[ii][3]

    ZP = ZP[0]

    print type(ZP), type(pfadnr), type(enigma), type(offset)
    print ZP, pfadnr, enigma, offset
    #ZP = '20180125170330'; pfadnr=0; enigma='neu'; offset=1


    # Pfade zu den Dateien festlegen
    #-----------------------------#

    year = ZP[0:4]; ye = ZP[2:4]; m = ZP[4:6]; d = ZP[6:8]; ht = ZP[8:10]; mt = ZP[10:12]; st = ZP[12:14]

    print('/automount/ags/velibor/gpmdata/dprV7/2A.GPM.DPR.V7-20170308.' + year + m + d + '*.HDF5')
    #pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.' + year + m + d + '*.HDF5')
    pfad_radar = glob.glob('/automount/ags/velibor/gpmdata/dprV7/2A.GPM.DPR.V7-20170308.' + year + m + d + '*.HDF5')

    print ('GPM:',pfad_radar)
    pfad_radar = pfad_radar[pfadnr]
    #pfad_radar_Ku = pfad_radar[0]

    deg_scan =  ["/ppi_1p5deg/","/ppi_2p4deg/","/ppi_3p4deg/",
                 "/n_ppi_010deg/","/n_ppi_045deg/",
                 "/n_ppi_082deg/","/n_ppi_110deg/","/n_ppi_140deg/",
                 "/n_ppi_180deg/","/n_ppi_280deg/","/n_vertical_scan/"][0]


    if enigma=='neu':
        print ('New enigma')
        deg_scan =  ["/n_ppi_010deg/"][0]

        boxpolpath = '/automount/radar/scans/' + year+ "/" +year +"-"+ m + "/" + year+ "-" + m +"-"+ d +\
                                   deg_scan+"*"+year+m+d+ht+mt+st+"*.h5"
        print (boxpolpath)
        ppi_datapath=glob.glob(boxpolpath)

        print ('Boxpol: ',ppi_datapath)
        ppi_datapath = ppi_datapath[0]


    else:
        try:
            ppi_datapath=glob.glob('/automount/radar-archiv/scans/' + year+ "/" +
                                   year +"-"+ m + "/" + year+ "-" + m +"-"+ d +
                                   deg_scan+ year + "-" + m +"-"+ d + "--" +ht +
                                   ":"+mt+":"+st+",*.mvol")
            print ppi_datapath
            ppi_datapath = ppi_datapath[0]

        except:
            ppi_datapath=glob.glob('/automount/radar/scans/' + year+ "/" +
                                   year +"-"+ m + "/" + year+ "-" + m +"-"+
                                   d + deg_scan+ year + "-" + m +"-"+ d +
                                   "--" +ht +":"+mt+":"+st+",*.mvol")
            print ('Old enigma')
            print ('Boxpol: ',ppi_datapath)
            ppi_datapath = ppi_datapath[0]


    # Wichtige Parameter festlegen
    TH = 15 #Threshold um Nullen fuer Niederschlag raus zu filtern


    #################################################### PPI BoxPol Daten einlesen
    #------------------------------------------------------------------------------

    ppi=h5py.File(ppi_datapath,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

    ZH0 = data['SCAN0']['ZH']['data']
    PHIDP = data['SCAN0']['PHIDP']['data']
    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']
    elevation=attrs['SCAN0']['elevation']
    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']
    rho = data['SCAN0']['RHOHV']['data']

    R = ZH0

    print ("________Beam Blockage______")
    R[151:165]=np.nan
    print (np.nanmax(R))

    print ("________CLUTTER______")
    rho_th  = 0.85
    R[rho<= rho_th] = np.nan
    print (np.nanmax(R))

    print ("________offset______")
    R = R + offset
    print (np.nanmax(R))
    #?
    print ("________ATTCORR______")

    """pia = wrl.atten.correctAttenuationHB(
        R,
        coefficients = dict(a=4.57e-5, b=0.731, gate_length=1.0),
        mode="warn",
        thrs=59.)
    pia[pia > 4.8] = 4.8

    print ("________ATTCORR2______")
    R = R + pia
    print (np.nanmax([R,pia]))"""

    print ("________DPR Threshold______")
    Z_boxpol = R

    ### Threshold for DPR sensitivity
    Z_boxpol[Z_boxpol<TH]=np.nan


    ################################################################# DPR Einlesen
    # -----------------------------------------------------------------------------
    scan = 'MS'
    gpmku = h5py.File(pfad_radar, 'r')
    gpmku_HS = gpmku[scan]['SLV']
    dpr_lat_1 = np.array(gpmku[scan]['Latitude'])
    dpr_lon_1 = np.array(gpmku[scan]['Longitude'])
    Z_dpr = np.array(gpmku_HS['zFactorCorrectedNearSurface'])
    Z_dpr[Z_dpr < TH] = np.nan

    ## Einlesen von Phase un Raintype
    dpr_raintype = np.array(gpmku[scan]['CSF']['typePrecip'], dtype=float)
    dpr_phase = np.array(gpmku_HS['phaseNearSurface'], dtype=float)


    ############################################################## RADOLAN einlesen
    # -----------------------------------------------------------------------------
    mtt = mt
    mtt = str(int(round(float(mtt)/5.0)*5.0))

    if mtt == '0':
        mtt = '00'
    if mtt == '5':
        mtt = '05'
    if mtt == '60':
        mtt = '55'

    r_pro = 'rx'

    pfad = ('/automount/radar/dwd/'+ r_pro +'/'+str(year)+'/'+str(year)+'-'+
            str(m)+'/'+ str(year)+'-'+str(m)+'-'+str(d)+'/raa01-'+r_pro+'_10000-'+
            str(ye)+str(m)+ str(d)+str(ht)+str(mtt)+'-dwd---bin.gz')

    print ('RADOLANPFAD: ', pfad)

    pfad_radolan = pfad[:-3]

    try:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad)
    except EnvironmentError:
        rw_filename = wradlib.util.get_wradlib_data_file(pfad_radolan)

    rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)

    radolan_zeit = rwattrs['datetime'].strftime("%Y.%m.%d -- %H:%M:%S")

    radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
    x = radolan_grid_xy[:,:,0]
    y = radolan_grid_xy[:,:,1]
    rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5

    ### Threshold for DPR sensitivity
    rwdata[rwdata<TH]=np.nan
    Z_radolan = rwdata


    ######################################################## Cut the Swath for Bonn
    # -----------------------------------------------------------------------------

    dpr_lon, dpr_lat, Z_dpr = cut_the_swath(dpr_lon_1,dpr_lat_1,Z_dpr, eu=0)
    dpr_lon, dpr_lat, dpr_raintype = cut_the_swath(dpr_lon_1,dpr_lat_1,dpr_raintype, eu=0)
    dpr_lon, dpr_lat, dpr_phase = cut_the_swath(dpr_lon_1,dpr_lat_1,dpr_phase, eu=0)

    ######################################################## Koordinaten Projektion
    # -----------------------------------------------------------------------------

    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    dpr_lon, dpr_lat = wradlib.georef.reproject(dpr_lon, dpr_lat,
                                                projection_target=proj_stereo ,
                                                projection_source=proj_wgs)

    blon, blat = wradlib.georef.reproject(blon0, blat0,
                                          projection_target=proj_stereo ,
                                          projection_source=proj_wgs)


    ############################################################### Dpr zuschneiden
    #------------------------------------------------------------------------------
    print ('Max --------------------------------------->', r[-1]/1000.)
    inner_r = 15. # in km!!!!

    print ('Min --------------------------------------->', inner_r)

    lon0, lat0, radius = blon, blat, r[-1]/1000.
    rr = np.sqrt((dpr_lat - lat0)**2 + (dpr_lon - lon0)**2)
    position = rr < radius

    # Maximum Range
    Z_dpr[np.where(rr > radius)] = np.nan
    dpr_raintype[np.where(rr > radius)] = np.nan
    dpr_phase[np.where(rr > radius)] = np.nan

    # Minimum Range
    Z_dpr[np.where(rr < inner_r)] = np.nan
    dpr_raintype[np.where(rr < inner_r)] = np.nan
    dpr_phase[np.where(rr < inner_r)] = np.nan


    ########################################################### RADOLAN zuschneiden
    #------------------------------------------------------------------------------

    rr2 = np.sqrt((y - lat0)**2 + (x - lon0)**2)
    position2 = rr2 < radius

    # Maximum Range
    Z_radolan[np.where(rr2 > radius)] = np.nan

    # Minimum Range
    Z_radolan[np.where(rr2 < inner_r)] = np.nan


    ########################################################### BoXPol zuschneiden
    #------------------------------------------------------------------------------
    # Minimum Range
    a = (inner_r*1000.)/(r[1]-r[0])# inner_r in m!!!!
    a = a.astype(int) # bin bis zu welchem radius alles NaN gesetzt wird
    Z_boxpol[:,0:a]=np.nan
    print (a)


    ############################################################### Inverse Dezibel
    #------------------------------------------------------------------------------
    Z_boxpol = wradlib.trafo.idecibel(Z_boxpol)
    Z_dpr = wradlib.trafo.idecibel(Z_dpr)
    Z_radolan = wradlib.trafo.idecibel(Z_radolan)

    ################################################## BoXPol interpolieren auf DPR
    #------------------------------------------------------------------------------
    t1 = clock()
    radar_location = (lon_ppi, lat_ppi, alt_ppi)
    #elevation = 1.5
    azimuths = az
    ranges = r
    polargrid = np.meshgrid(ranges, azimuths)
    lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1],
                                                     elevation, radar_location)
    lon, lat = wradlib.georef.reproject(lon, lat, projection_target=proj_stereo ,
                                        projection_source=proj_wgs)

    grid_xy = np.vstack((dpr_lon.ravel(), dpr_lat.ravel())).transpose()

    xy=np.concatenate([lon.ravel()[:,None],lat.ravel()[:,None]], axis=1)



    Z_boxpol_ipoli, Z_boxpol_ipoli_std, Z_boxpol_ipoli_median, Z_boxpol_ipoli_max, Z_boxpol_ipoli_min = ipoli_radi_stat(xy,Z_boxpol.ravel(),grid_xy,2.5)

    Z_boxpol_ipoli = Z_boxpol_ipoli.reshape(dpr_lon.shape)
    Z_boxpol_ipoli_std = Z_boxpol_ipoli_std.reshape(dpr_lon.shape)
    Z_boxpol_ipoli_median = Z_boxpol_ipoli_median.reshape(dpr_lon.shape)
    Z_boxpol_ipoli_max = Z_boxpol_ipoli_max.reshape(dpr_lon.shape)
    Z_boxpol_ipoli_min = Z_boxpol_ipoli_min.reshape(dpr_lon.shape)

    t2 = clock()
    print ('Interpolationsdauer Boxpol auf DPR:', t2 - t1)


    ################################################# RADOLAN interpolieren auf DPR
    #------------------------------------------------------------------------------

    xy_rado = np.vstack((x.ravel(), y.ravel())).transpose()

    Z_radolan_ipoli, Z_radolan_ipoli_std, Z_radolan_ipoli_median, Z_radolan_ipoli_max, Z_radolan_ipoli_min = ipoli_radi_stat(xy_rado,Z_radolan.ravel().filled(np.nan),grid_xy,2.5)

    Z_radolan_ipoli = Z_radolan_ipoli.reshape(dpr_lon.shape)
    Z_radolan_ipoli_std = Z_radolan_ipoli_std.reshape(dpr_lon.shape)
    Z_radolan_ipoli_median = Z_radolan_ipoli_median.reshape(dpr_lon.shape)
    Z_radolan_ipoli_max = Z_radolan_ipoli_max.reshape(dpr_lon.shape)
    Z_radolan_ipoli_min = Z_radolan_ipoli_min.reshape(dpr_lon.shape)

    t3 = clock()
    print ('Interpolationsdauer RADOLAN auf DPR:', t3 - t2)

    ################################################################### in  Dezibel
    #------------------------------------------------------------------------------
    Z_boxpol = wradlib.trafo.decibel(Z_boxpol)
    Z_dpr = wradlib.trafo.decibel(Z_dpr)
    Z_radolan = wradlib.trafo.decibel(Z_radolan)
    Z_boxpol_ipoli = wradlib.trafo.decibel(Z_boxpol_ipoli)
    Z_radolan_ipoli = wradlib.trafo.decibel(Z_radolan_ipoli)



    #zwischenspeichern
    #np.save('/automount/ags/velibor/data/dpr_boxpol_radolan/test/median/'+ZP+'.npy',
    #        [Z_dpr, Z_boxpol, Z_boxpol_ipoli, Z_radolan, Z_radolan_ipoli, dpr_phase, dpr_raintype, dpr_lon, dpr_lat])

tend = clock()
print ('Gesamtlaufzeit: ',(tend-tstart)/60., 'min')

