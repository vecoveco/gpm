 
'''Dieses Program soll dazu dienen die 
Radardaten von BoxPol mit den GPM Daten
hinsichtlich der Regenraten zu validieren.
Hier mehrere Ueberfluge und mehrere
 interpolations Methoden analysiert'''

#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------
'''Einlesen von Modulen'''
#--------------------------------------------------------------------------------------------------------
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob

'''ftp://ftp.meteo.uni-bonn.de/pub/pablosaa/gpmdata/'''

#--------------------------------------------------------------------------------------------------------
### Pfad mit String ##
#--------------------------------------------------------------------------------------------------------
TH = 0.01 #Threshold um Nullen fuer Niederschlag raus zu filtern

C_Idw=[]

C_Linear=[]

C_Nearest=[]

C_OrdKriging=[]

C_Idw_std=[]

C_Linear_std=[]

C_Nearest_std=[]

C_OrdKriging_std=[]

ipoli = [wradlib.ipol.Idw, wradlib.ipol.Linear, wradlib.ipol.Nearest, wradlib.ipol.OrdinaryKriging]

#LISTE der Ueberfluege des GPM mit Niederschlagevents
#BITTE ZEITRAUM DES PPI EINTRAGEN!
LISTE = ("20141007023000","20150128172208","20140629145000","20140629145925","20140826145000","20140826150322","20140921070500","20140921071058","20141007023744","20141016001500","20141016002458","20150128171500")#
LISTE=sorted(LISTE)

for i in range(0,len(LISTE)):
    for j in range(0,len(ipoli)):  
        ZP=LISTE[i]
        print "initialize ::: " + ZP
    
        year = ZP[0:4]
        m = ZP[4:6]
        d = ZP[6:8]
        ht = ZP[8:10]
        mt = ZP[10:12]
        st = ZP[12:14]



        pfad = ('/home/velibor/shkgpm/data/' + year + m + d + '/radar/*.HDF5')
        pfad_radar = sorted(glob.glob(pfad))
        pfad_radar_Ku = pfad_radar[2]
        pfad_radar_DPR = pfad_radar[0]
        pfad_radar_DPRGMI = pfad_radar[3]
        pfad_radar_Ka = pfad_radar[1]

        pfad1 = ('/home/velibor/shkgpm/data/' + year + m + d + '/boxpol/gpm_rhi_01/*.mvol')
        pfad_boxpol = glob.glob(pfad1)
        pfad_boxpol_rhi01 = pfad_boxpol[0]

        pfad2 = ('/home/velibor/shkgpm/data/' + year + m + d + '/gprof/*.HDF5')
        pfad_gprof = glob.glob(pfad2)
        pfad_gprof_g = pfad_gprof[0]

        ppi_datapath=('/automount/radar/scans/' + year+ "/" + year +"-"+ m + "/" + year+ "-" + m +"-"+ d + "/ppi_1p5deg/"+ year + "-" + m +"-"+ d + "--" +ht +":"+mt+":"+st+",00.mvol")

#--------------------------------------------------------------------------------------------------------
### Rhi BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

        gpmrhi01 = h5py.File(pfad_boxpol_rhi01, 'r')

#--------------------------------------------------------------------------------------------------------
### PPI BoxPol Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

        ppi=h5py.File(ppi_datapath,'r')
        data, attrs = wradlib.io.read_GAMIC_hdf5(ppi_datapath)

        ZH = data['SCAN0']['ZH']['data']
        r = attrs['SCAN0']['r']
        az = attrs['SCAN0']['az']
        lon_ppi = attrs['VOL']['Longitude']
        lat_ppi = attrs['VOL']['Latitude']
        alt_ppi = attrs['VOL']['Height']
#Umwandeln von Z in RR
        Z = wradlib.trafo.idecibel(ZH)
        R = wradlib.zr.z2r(Z, a=200., b=1.6) # Rainfall intensity using the Marshal-Palmer Z(R) parameters
        depth = wradlib.trafo.r2depth(R, 300)# Convert to rainfall depth assuming a duration of 300 seconds
#--------------------------------------------------------------------------------------------------------
### gprof Daten einlesen ##
#--------------------------------------------------------------------------------------------------------

        gpmgmi = h5py.File(pfad_gprof_g, 'r')


#--------------------------------------------------------------------------------------------------------
###---- GPROF ---- Einlesen von Oberflachenniederschlag und Lat/Lon von Gprof ##
#--------------------------------------------------------------------------------------------------------

        gpmgmi.keys()
        gpmgmi_S1=gpmgmi['S1']				
        gprof_lat=gpmgmi_S1['Latitude']			#(2962, 221)
        gprof_lon=gpmgmi_S1['Longitude']			#(2962, 221)
        gprof_pp=gpmgmi_S1['surfacePrecipitation']   	#(2962, 221)


#--------------------------------------------------------------------------------------------------------
### ------- In Arrays umwandeln ----------- ##
#--------------------------------------------------------------------------------------------------------

        gprof_pp_a = np.array(gprof_pp)
        gprof_lon_a = np.array(gprof_lon)
        gprof_lat_a = np.array(gprof_lat)

#--------------------------------------------------------------------------------------------------------
### ------- Lon Lat Bestimmung ----------- ##
#--------------------------------------------------------------------------------------------------------
# Indexzugriff
#   gprof_lat_a[[Zeile][Zeilenelement]]
# Indexsuche
#   itemin= np.where((gprof_lat_a<58) & (gprof_lat_a>57))
        bonn_lat1 = 49.9400
        bonn_lat2 = 51.3500
        bonn_lon1 = 6.40000
        bonn_lon2 = 8.10000

        ilat= np.where((gprof_lat_a>49.9400) & (gprof_lat_a<51.3500))
        ilon= np.where((gprof_lon_a>6.40000) & (gprof_lon_a<8.10000))
        lonstart = ilon[0][0]	#erstes Element
        lonend = ilon[0][-1]	#letztes Element
        latstart = ilat[0][0]
        latend = ilat[0][-1]



        gprof_pp_a[gprof_pp_a==-9999] = np.nan #NAN!
        radar_location = (lon_ppi, lat_ppi, alt_ppi) # (lon, lat, alt) in decimal degree and meters
        elevation = 1.5 # in degree
        azimuths = az # in degrees
        ranges = r # in meters
        polargrid = np.meshgrid(ranges, azimuths)
        lon, lat, alt = wradlib.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation, radar_location)

        gk3 = wradlib.georef.epsg_to_osr(31467)
        x, y = wradlib.georef.reproject(lon, lat, projection_target=gk3)
        xgrid, ygrid = wradlib.georef.reproject(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], projection_target=gk3)

        grid_xy = np.vstack((xgrid.ravel(), ygrid.ravel())).transpose()

        xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
    
        
	gridded = wradlib.comp.togrid(xy, grid_xy, ranges[-1], np.array([x.mean(), y.mean()]), R.ravel(), ipoli[j])#Linear, Idw, Nearest,OrdinaryKriging,
	gridded = np.ma.masked_invalid(gridded).reshape(xgrid.shape)
    #PLOT
        fig = plt.figure(figsize=(13,10))
#Levels berechnen
	maxv = np.max([np.max(gridded),np.max(np.ma.masked_invalid(gprof_pp_a)[latstart:latend])])
#
	plt.subplot(221)

	A = gridded
	B = np.ma.masked_invalid(gprof_pp_a)[latstart:latend]
#Nullen entfernen
	A[A<TH]=np.nan
	B[B<TH]=np.nan


#Scatter mit regrssion
	from scipy import stats
	mask = ~np.isnan(B) & ~np.isnan(A)
	slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
	line = slope*B+intercept
	plt.plot(B,line,'r-',B,A,'o')
#plt.gca().set_aspect('equal', adjustable='box')		#control quadratic x y axis
#plt.axis('equal')						#control quadratic x y axis
	plt.xlim(0,maxv)						#control quadratic x y axis
	plt.ylim(0,maxv)						#control quadratic x y axis
	plt.xlabel("gprof")
	plt.ylabel("ppi BoxPol")
	plt.title("Scatterplot (gprof/ppi), cor: " + str(r_value))


##
	plt.subplot(222)
#rainrate
	ax1, pm2 = wradlib.vis.plot_ppi(R,r,az,vmin=0,vmax=maxv)
	cbar = plt.colorbar(pm2, shrink=0.75)
	cbar.set_label("RR [mm/h]")
	plt.xlim((-101000,101000))
        plt.ylim((-101000,101000))
	plt.xticks(())
	plt.yticks(())
	plt.xlabel("X Range [km]")
	plt.ylabel("Y Range [km]")
	plt.title(ppi_datapath[-28:-8])

	plt.subplot(223)
#gprof rainrate
	pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], np.ma.masked_invalid(gprof_pp_a)[latstart:latend],vmin=0,vmax=maxv)
	plt.xlim((bonn_lon1,bonn_lon2))
	plt.ylim((bonn_lat1,bonn_lat2))
	plt.title(pfad_boxpol_rhi01[-28:-6])
	cbar = plt.colorbar(pm2, shrink=.75)
	cbar.set_label("GPROF surfacePrecipitation [mm/h]")
	plt.xlabel("Easting (m)")
	plt.ylabel("Northing (m)")

	plt.subplot(224)
#ppi rainrate
	pm2 = plt.pcolormesh(gprof_lon_a[latstart:latend], gprof_lat_a[latstart:latend], gridded,vmin=0,vmax=maxv)#xgrid, ygrid,
#plt.xlim((lon.min(),lon.max()))
#plt.ylim((lat.min(),lat.max()))
	plt.xlim((bonn_lon1,bonn_lon2))
	plt.ylim((bonn_lat1,bonn_lat2))
	plt.title(ppi_datapath[-28:-8])
	cbar = plt.colorbar(pm2, shrink=0.75)
	cbar.set_label("Boxpol_ppi_interpolation RR [mm/h]")
	plt.xlabel("Easting (m)")
	plt.ylabel("Northing (m)")

	plt.savefig('plot/' + ppi_datapath[-28:-8] + '_Vergleich_ohneNull_' +str(ipoli[j]) +'.png')
	plt.close()
    
    
    
	if j == 0:
            C_Idw.append(r_value)
            C_Idw_std.append(std_err)
        elif j == 1:    
            C_Linear.append(r_value)
            C_Linear_std.append(std_err)
        elif j == 2:
            C_Nearest.append(r_value)
            C_Nearest_std.append(std_err)
        elif j == 3:
            C_OrdKriging.append(r_value)
            C_OrdKriging_std.append(std_err)



df = pd.DataFrame({'Idw':C_Idw,'Linear':C_Linear,'Nearest':C_Nearest,'OrdKriging':C_OrdKriging,'Idw_std':C_Idw_std,'Linear_std':C_Linear_std,'Nearest_std':C_Nearest_std,'OrdKriging_std':C_OrdKriging_std})
ts = pd.to_datetime(LISTE)
df=df.set_index(ts) #Zeitrheie erstellen und einfuegen

df[['Idw','Linear','Nearest','OrdKriging']].plot(kind="bar")
plt.xticks(fontsize=8,rotation=23)
plt.title("Correlation boxpol/gprof [All_ipol]")
plt.ylabel("Correlation")
plt.grid()
plt.savefig("plot/BAR_Corr.png")
plt.close()

#--------------------------------------------------------------------------------------------------------
### ------- Alle Korrelationen ----------- ##
#--------------------------------------------------------------------------------------------------------
'''
df_corr= df[['Idw','Linear','Nearest','OrdKriging']]
df_std = df[['Idw_std','Linear_std','Nearest_std','OrdKriging_std']]
fig, ax = plt.subplots()
df_corr.plot(yerr=df_std, ax=ax, kind='box')
plt.savefig("plot/00PLOT.png")
'''
'''
y_LISTE = np.arange(len(LISTE))
plt.barh(y_LISTE,Corre,align='center', alpha=0.4)
plt.yticks(y_LISTE, LISTE,fontsize=7)
plt.xlim(-1, 1)
plt.xlabel('Correlation')
plt.title('Correlation_boxpol_gprof')
plt.grid()
plt.savefig("CORRE_ALL.png")

plt.close()
'''
