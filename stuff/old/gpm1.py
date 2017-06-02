'''Dieses Program soll dazu dienen die 
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Regenraten zu validieren.
Hier wird nur ein Ueberflug analysiert'''

#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------
'''Einlesen von Modulen'''
#--------------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import cKDTree

def mask_from_bbox(x, y, bbox):
    """Return 2-d index array based on spatial selection from a bounding box.

    Use this function to create a 2-d boolean mask from 2-d arrays of grids points.

    Parameters
    ----------
    x : nd array of shape (num rows, num columns)
        x (Cartesian) coordinates
    y : nd array of shape (num rows, num columns)
        y (Cartesian) coordinates
    bbox : dictionary with keys "left", "right", "bottom", "top"
        These must refer to the same Cartesian reference system as x and y

    Returns
    -------
    out : mask, shape
          mask is a boolean array that is True if the point is inside the bbox
          shape is the shape of the True subgrid

    """
    ny, nx = x.shape

    ix = np.arange(x.size).reshape(x.shape)

    # Find bbox corners
    #    Plant a tree
    tree = cKDTree(np.vstack((x.ravel(),y.ravel())).transpose())
    # find lower left corner index
    dists, ixll = tree.query([bbox["left"], bbox["bottom"]], k=1)
    ill = (ixll / nx)-1
    jll = (ixll % nx)-1
    # find upper right corner index
    dists, ixur = tree.query([bbox["right"],bbox["top"]], k=1)
    iur = (ixur / nx)+1
    jur = (ixur % nx)+1

    mask = np.repeat(False, ix.size).reshape(ix.shape)

    if iur>ill:
        mask[ill:iur,jll:jur] = True
        shape = (iur-ill, jur-jll)
    else:
        mask[iur:ill,jll:jur] = True
        shape = (ill-iur, jur-jll)
    return mask, shape


#Einlesen mit String 
#probe="/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_01/2014-06-29--14:54:52,00.mvol"
#probe = h5py.File(probe, 'r')

'''ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/'''
#--------------------------------------------------------------------------------------------------------
### Pfad mit String ##
#--------------------------------------------------------------------------------------------------------
'''
y = '2014'
m = '06'
d = '29'
ht = '14'
mt = '54'
st = '52'
Skoo = '140300'
Ekoo = '153533.001897'

pfad_rhi = '/home/velibor/shkgpm/data/' + y + m + d + '/boxpol/gpm_rhi_01/' + y + '-' + m + '-' + d + '--' + ht + ':' + mt + ':' + st + ',00.mvol'
pfad_gmi='/home/velibor/shkgpm/data/' + y + m + d + '/gprof/2A.GPM.GMI.GPROF' + y + 'v1-4.' + y + m + d + '-S'+ Skoo +'-E'+ Ekoo +'.V03C.HDF5'
'''
#--------------------------------------------------------------------------------------------------------
### Hier kann man gucken welchen Pfadtyp man dann zum einlesen benutzt, den oben oder den unten ##
#--------------------------------------------------------------------------------------------------------

gpmrhi01_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_01/2014-06-29--14:54:52,00.mvol'
gpmrhi02_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_02/2014-06-29--14:54:52,00.mvol'
gpmrhi03_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_03/2014-06-29--14:54:52,00.mvol'
gpmrhi04_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_04/2014-06-29--14:54:52,00.mvol'
gpmrhi05_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_05/2014-06-29--14:54:52,00.mvol'
gpmrhi06_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_06/2014-06-29--14:54:52,00.mvol'
gpmrhi07_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_07/2014-06-29--14:54:52,00.mvol'
gpmrhi08_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_08/2014-06-29--14:54:52,00.mvol'
gpmrhi09_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_09/2014-06-29--14:54:52,00.mvol'
gpmrhi10_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_10/2014-06-29--14:54:52,00.mvol'
gpmrhi11_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_11/2014-06-29--14:54:52,00.mvol'
gpmrhi12_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_12/2014-06-29--14:54:52,00.mvol'
gpmrhi13_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_13/2014-06-29--14:54:52,00.mvol'
gpmrhi14_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_14/2014-06-29--14:54:52,00.mvol'
gpmrhi15_datapath = '/home/velibor/shkgpm/data/20140629/boxpol/gpm_rhi_15/2014-06-29--14:54:52,00.mvol'

gpmgmi_datapath = '/home/velibor/shkgpm/data/20140629/gprof/2A.GPM.GMI.GPROF2014v1-4.20140629-S140300-E153533.001897.V03C.HDF5'
dpr_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.DPR.V5-20140827.20140629-S140300-E153533.001897.V03B.HDF5'
ka_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.Ka.V5-20140829.20140629-S140300-E153533.001897.V03B.HDF5'
ku_datapath = '/home/velibor/shkgpm/data/20140629/radar/2A.GPM.Ku.V5-20140829.20140629-S140300-E153533.001897.V03B.HDF5'
dprgmi_datapath = '/home/velibor/shkgpm/data/20140629/radar/2B.GPM.DPRGMI.CORRA2014.20140629-S140300-E153533.001897.V03C.HDF5'

#--------------------------------------------------------------------------------------------------------
### BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

gpmrhi01 = h5py.File(gpmrhi01_datapath, 'r')
gpmrhi02 = h5py.File(gpmrhi02_datapath, 'r')
gpmrhi03 = h5py.File(gpmrhi03_datapath, 'r')
gpmrhi04 = h5py.File(gpmrhi04_datapath, 'r')
gpmrhi05 = h5py.File(gpmrhi05_datapath, 'r')
gpmrhi06 = h5py.File(gpmrhi06_datapath, 'r')
gpmrhi07 = h5py.File(gpmrhi07_datapath, 'r')
gpmrhi08 = h5py.File(gpmrhi08_datapath, 'r')
gpmrhi09 = h5py.File(gpmrhi09_datapath, 'r')
gpmrhi10 = h5py.File(gpmrhi10_datapath, 'r')
gpmrhi11 = h5py.File(gpmrhi11_datapath, 'r')
gpmrhi12 = h5py.File(gpmrhi12_datapath, 'r')
gpmrhi13 = h5py.File(gpmrhi13_datapath, 'r')
gpmrhi14 = h5py.File(gpmrhi14_datapath, 'r')
gpmrhi15 = h5py.File(gpmrhi15_datapath, 'r')

#--------------------------------------------------------------------------------------------------------
### gprof Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

gpmgmi = h5py.File(gpmgmi_datapath, 'r')

#--------------------------------------------------------------------------------------------------------
### radar Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

dpr = h5py.File(dpr_datapath, 'r')
ka = h5py.File(ka_datapath, 'r')
ku = h5py.File(ku_datapath, 'r')
dprgmi = h5py.File(dprgmi_datapath, 'r')

#--------------------------------------------------------------------------------------------------------
###---- GPROF ---- Einlesen von Oberflachenniederschlag und Lat/Lon von Gprof ##
#--------------------------------------------------------------------------------------------------------

gpmgmi.keys()
gpmgmi_S1=gpmgmi['S1']				
gprof_lat=gpmgmi_S1['Latitude']			#(2962, 221)
gprof_lon=gpmgmi_S1['Longitude']		#(2962, 221)
gprof_pp=gpmgmi_S1['surfacePrecipitation']   	#(2962, 221)

#--------------------------------------------------------------------------------------------------------
### ---- Rhi BoXPol ---- Einlesen... ##
#--------------------------------------------------------------------------------------------------------

scan=gpmrhi01['scan0']
boxpol_zdr=scan['moment_9']   			#(58, 1800)


#--------------------------------------------------------------------------------------------------------
### ---- radar dpr ka ku dprgmi---- ##
#--------------------------------------------------------------------------------------------------------

dpr_NS = dpr['NS']				#Einlesen ueber HS und MS
dpr_SLV = dpr_NS['SLV']
dpr_pr = dpr_SLV['precipRate']			#(7934, 49, 176)
dpr_prs = dpr_SLV['precipRateNearSurface']	#(7934, 49)
dpr_lat = dpr_NS['Latitude']			#(7934, 49)
dpr_lon = dpr_NS['Longitude']			#(7934, 49)


#--------------------------------------------------------------------------------------------------------
### ------- In Arrays umwandeln ----------- ##
#--------------------------------------------------------------------------------------------------------

gprof_pp_a = np.array(gprof_pp)
gprof_lon_a = np.array(gprof_lon)
gprof_lat_a = np.array(gprof_lat)

dpr_lon_a = np.array(dpr_lon)
dpr_lat_a = np.array(dpr_lat)
dpr_pr_a = np.array(dpr_pr)
dpr_prs_a = np.array(dpr_prs)



#--------------------------------------------------------------------------------------------------------
### ------- Befehl ----------- ##
#--------------------------------------------------------------------------------------------------------
# Indexzugriff
#   gprof_lat_a[[Zeile][Zeilenelement]]
# Indexsuche
#   itemin= np.where((gprof_lat_a<58) & (gprof_lat_a>57))

bbox = dict(left=6.4, right=8.1, bottom=49.94, top=51.35)

index=np.argsort(gprof_lon_a.ravel(),axis=0)
gprof_lon_a_sort = gprof_lon_a.ravel()[index]
gprof_lat_a_sort = gprof_lat_a.ravel()[index]  
gprof_pp_a_sort = gprof_pp_a.ravel()[index] 
#mask, shape = mask_from_bbox(gprof_lon_a_sort, gprof_lat_a_sort, bbox)
ilon= np.where((gprof_lon_a_sort>6.40000) & (gprof_lon_a_sort<8.10000))

ilat= np.where((gprof_lat_a_sort>49.9400) & (gprof_lat_a_sort<51.3500))

in1d = np.where(np.in1d(ilon, ilat))[0]

#mask = ilon[in1d]


#lonstart = ilon[0][0]#erstes Element
#lonend = ilon[0][-1]#letztes Element
#latstart = ilat[0][0]
#latend = ilat[0][-1]
#--------------------------------------------------------------------------------------------------------
''' --------------------------------- PLOTS ------------------------------------- ''' 
#--------------------------------------------------------------------------------------------------------


### Plot GPROF pcolormesh
#plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], gprof_pp_a[latstart:latend],clim=(0.0, 70))
plt.pcolormesh(gprof_lon_a_sort[mask], gprof_lat_a_sort[mask], gprof_pp_a_sort[mask])#,clim=(0.0, 70)
plt.savefig('gprof_lonlat_data.png')
plt.close()
#

#yAchse eingrenzen auf Bonn
plt.pcolormesh(gprof_lon_a_sort, gprof_lat_a_sort, gprof_pp_a_sort)#,clim=(0.0, 70)
#plt.ylim((25,250))
plt.savefig('gprof_lonlat_data.png')
plt.close()


plt.imshow(gprof_pp, interpolation='nearest', cmap='Blues', origin='lower',clim=(0.0, 2)) #clim veraendert colorbar!!!!!!
plt.colorbar(shrink=.92)
#plt.xticks(())		#entfernt x ticks
#plt.yticks(())		#entfernt y ticks
plt.savefig('gprof_pp1.png')
plt.close()


plt.imshow(boxpol_zdr, interpolation='nearest', cmap='Blues', origin='lower')
plt.colorbar(shrink=.92)
plt.savefig('boxpol_pp.png')
plt.close()




plt.imshow(dpr_prs, interpolation='nearest', cmap='Blues', origin='lower')
plt.colorbar(shrink=.92)
plt.savefig('dpr_prs.png')
plt.close()

plt.imshow(dpr_lat, interpolation='nearest', cmap='Blues', origin='lower')
plt.colorbar(shrink=.92)
plt.savefig('dpr_lat.png')
plt.close()

#plt.plot( dpr_lat_a, dpr_lon_a,dpr_pr_a)
#plt.colorbar()
#plt.savefig('TEiST.png')
#plt.close()
