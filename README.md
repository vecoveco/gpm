# Global Precipitation Measuring Validation
============================================
![pic](https://pmm.nasa.gov/sites/default/files/imce/GPM_banner_3_droplets.jpg) 
> The Global Precipitation Measurement (GPM) mission is an international network of satellites that provide the next-generation global observations of rain and snow. 
> [pmm.nasa.gov]


## Description:
----------------

In order to understand the global distribution of precipitation as widely as possible, station- and radar measurements 
are not sufficient. Considering the meteorological observation network clearly shows that a direct observation on water
surfaces is virtually absent on the basis of conventional station equipment. In polar and desert regions this is also true.
Satellites could solve this problem. For this purpose, the GPM-satellite is equipped with a microwave camera
(GMI GPM Microwave Imager) and a dual-frequency radar (DPR Dual-frequency Precipitation Radar). 
Following the long collaboration participation of the University of Bonn at several NASA GPM ground-validation field-campaigns,
the polarimetric research radar BoXPol has became one of the few European research radars part of an international network
for ground validation of the GPM satellite. Routinely BoXPol performs specific scans for GPM overpasses, and its polarimetric
capabilities are exploited in order to estimate high quality rain rate fields within the GPM overpass field-of-view. 
Therefore the BoxPol retrievals are used to evaluate the corresponding products from the different suit of algorithms
currently used by GPM. For each overpass of GPM on the region detected by the Bonner radar, two and three dimensional data is created. This provides information about the current precipitation. To compare the GPM and BoXPol data, the polar data of the Bonn radar have to be interpolated on the GPM grid.

[References by NASA](http://svs.gsfc.nasa.gov/Gallery/GPM.html)

[References by JAXA](http://global.jaxa.jp/projects/sat/gpm/)



folgt...

## Application:
---------------
These scripts are used to compare satellites and radar data and provide a statement about their quality or correlation.

Verwendung folgt ...


## Scripts and 	program structure:
--------------------------------
folgt...


## Satellite Information:
------------------------------

Satellit | GPM Global Precipitation Measuring
--------------|---
Typ |  Erdoberflächensatellit
Manufacturer | NASA JAXA
Launch date | 27.02.2014
Launch site | Space Center Tanegashima Japan
Path height | 407 km
Orbital period | 92.6 min
Path inclination | 65 Grad
Contractor  |  Mitsubishi
Rocket  | 	H-IIA 202
Reference system |	Geocentric
Orbit Type |	Low Earth Orbit LEO

Instrument | Abk. | Spur | Typ  | Produkt
------------|-----|------|------|----------
Microwave Imager | GMI | 885 km | 13 channels | 2dGPROF
Dual_Frequency Precipitation Radar | DPR Ka-Band | 120 km | 35.5 GHz | 3dKa
Dual_Frequency Precipitation Radar | DPR Ku-Band | 245 km | 13.6 GHz | 3dKu



Satellite | Agency | Short | Nation
----------|---------|--------| -------
Megha-Tropiques | Centre National d’Études Spatiales,  Indian Space Research Organization  | CNES, ISRO | France, India
NOAA 18 |  National Oceanic and Atmospheric Administration | NOAA | USA
NOAA 19 |  National Oceanic and Atmospheric Administration | NOAA | USA
GCOM_W1 |  Japan Aerospace Exploration Agency | JAXA | Japan
DMSP F17 |  Department of Defense | DOD | USA
DMSP F18 |  Department of Defense | DOD | USA
DMSP F19 |  Department of Defense | DOD | USA
DMSP F20 |  Department of Defense | DOD | USA
JPSS-1 | National Oceanic and Atmospheric Administration | NOAA | USA
TRMM |  National Aeronautics and Space Administration ,Japan Aerospace Exploration Agency | NASA, JAXA | USA, Japan
MetOp B |  European Organization for the Exploitation of Meteorological Satellites | EUMETSAT | Europa
MetOp C |  European Organization for the Exploitation of Meteorological Satellites | EUMETSAT | Europa
Suomi NPP | National Aeronautics and Space Administration, National Oceanic and Atmospheric Administration | NASA, NOA | USA




==================
# GPM Public Data
==================

> ftp://arthurhou.pps.eosdis.nasa.gov/

/gpmdata contains the latest version of available data products. 

Directories are laid out as:

*/gpmdata/documents* !!!! Please review the contents of 
                       this directory. It contains important 
                       data notices and caveats !!!!
                       
                       
*/gpmdata/geolocation*  
> contains geolocation related files

*/gpmdata/coincidence* 
> contains satellite-ground coincidence information

*/gpmdata/browse* 
> contains PNG browse images of products.

Data products exist under:
*/gpmdata/YYYY/MM/DD/*

 > base/reg  - 1B base radiometer products
 
 > 1B    - 1B radiometer products for GPMcore and all partners
 
 > 1C    - 1C radiometer products for GPMcore and all partners
 
 > radar - L2 and L3 products from DPR and Combined
 
 > gprof - L2 and L3 products from GPROF for GPMcore and all partners
 
 > imerg - IMERG products
 
 
 
---------------------- 
# Level
----------------------
 
 **Level 1**
 > Level 1A: Reconstructed, unprocessed instrument data at full resolution, time referenced, and annotated with ancillary information, including radiometric and geometric calibration coefficients and georeferencing parameters (i.e., platform ephemeris), computed and appended, but not applied, to Level 0 data.

 > Level 1B: Radiometrically corrected and geolocated Level 1A data that have been processed to sensor units..

 > Level 1C: Common intercalibrated brightness temperature (Tc) products using the GPM Microwave Imager (GMI) Level 1B as the reference standard.
 
 **Level 2**
 
  > Derived geophysical parameters at the same resolution and location as those of the Level 1 data.
 
 **Level 3**
 
  > Geophysical parameters that have been spatially and/or temporally resampled from Level 1 or Level 2 data.
 
 
---------------------- 
# Products
---------------
 
 ----------------------------------------------------------
### 1A
----------------------------------------------------------

> 1A.GPM.GMI.COUNT2016.20150609-S000921-E014155.007256.V04A.HDF5

> 1A-GMI: GMI Packet Data Transmitted by the Satellite

> Level 1A: Reconstructed, unprocessed instrument data at full resolution, time referenced, and annotated with ancillary information, including radiometric and geometric calibration coefficients and georeferencing parameters (i.e., platform ephemeris), computed and appended, but not applied, to Level 0 data.

esolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
4km x 4km - 16 orbits per day| 	Latitudes 70°N-S, March 2014 - present| 	40 hours (Prod)| 	HDF5 |	Prod: FTP (PPS)*(/YYYY/MM/DD/1A)

----------------------------------------------------------
### 1B
----------------------------------------------------------

> 1B.GPM.GMI.TB2015.20150609-S000921-E014155.007256.V04A.HDF5

> 1B-GMI: GMI brightness temperatures

> The Level 1B algorithm and software transform Level 0 counts into geolocated and calibrated antenna temperatures (Ta) and brightness temperatures (Tb). Ta is obtained by utilizing the sensor radiometric calibration as well as various corrections based on after launch analyses. Tb is derived from Ta after antenna pattern correction (APC) and along scan corrections. Figure 1.16 describes the relationship between algorithm components and products (or output). 

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
Varies by Channel - 16 orbits per day |	Latitudes 70°N-S, Past 2 Weeks (NRT) |	20 minutes (NRT); 6 hours (Prod) 	|HDF5 |	FTP (PPS) (/YYYY/MM/DD/1B)


----------------------------------------------------------
### 1C
----------------------------------------------------------

> 1C.GPM.GMI.XCAL2015-C.20150609-S000921-E014155.007256.V04A.HDF5

> 1C-GMI: Calibrated GMI brightness temperatures

> The Level 1C algorithms transform equivalent Level 1B radiance data into Level 1C products. The input source data are geolocated and radiometric calibrated antenna temperature (Ta) or brightness temperature (Tb). The output Level 1C products are common intercalibrated brightness temperature (Tc) products using the GPM Microwave Imager (GMI) as the reference standard. The Level 1C algorithms contain the following major components:

    Orbitization.
    Satellite intercalibration.
    Quality control.
    Ancillary data calculations.

The detail of L1C algorithms and implementation depends on the details of each sensor. In this document, the Level 1C algorithms are described in a general sense. Individual sensor-specific details are provided separately in Appendices A through G: A) GMI, B) LIC-R GMI, C) Tropical Rainfall Measuring Mission (TRMM) Microwave Imager (TMI), D) Special Sensor Microwave Imager/Sounder (SSMI/S), E) Advanced Microwave Scanning Radiometer 2 (AMSR2), F) Advanced Technology Microwave Sounder (ATMS), G) Sondeur Atmospherique du Profil d'Humidite Intertropicale par Radiometrie (SAPHIR), and H) Microwave Humidity Sounder (MHS).

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
Varies by Channel - 16 orbits per day |	orbital, Past 2 Weeks (NRT) |	20 minutes (NRT); 6 hours (Prod) 	|HDF5 |	FTP (PPS) (/YYYY/MM/DD/1C)

----------------------------------------------------------
### gprof
----------------------------------------------------------

> 2A.GPM.GMI.GPROF2014v2-0.20150609-S000921-E014155.007256.V04A.HDF5

> 2A-GPROF-GMI: GMI single-orbit rainfall estimates

> The 2AGPROF (also known as, GPM GPROF (Level 2)) algorithm retrieves consistent precipitation and related science fields from the following GMI and partner passive microwave sensors: GMI, SSMI (DMSP F15), SSMIS (DMSP F16, F17, F18) AMSR2 (GCOM-W1), TMI MHS (NOAA 18&19, METOP A&B), ATMS (NPP), SAPHIR (MT1) This provides the bulk of the 3-hour coverage achieved by GPM. For each sensor, there are near-realtime (NRT) products, standard products, and climate products. These differ only in the amount of data that are available within 3 hours, 48 hours, and 3 months of collection, as well as the ancillary data used. 

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
4km x 4km, 16 orbits per day |	Latitudes 90°N-S, Past 2 weeks (NRT); March 2014 - present (Prod) |	2 hours (NRT); 40 hours (Prod)	|HDF5 |	FTP (PPS) (/YYYY/MM/DD/2A.GPM.GMI)

----------------------------------------------------------
### imerg
----------------------------------------------------------

> 3B-HHR.MS.MRG.3IMERG.20150609-S000000-E002959.0000.V03D.HDF5

> IMERG: Rainfall estimates combining data from all passive-microwave instruments in the GPM Constellation

> This algorithm is intended to intercalibrate, merge, and interpolate “all” satellite microwave precipitation estimates, together with microwave-calibrated infrared (IR) satellite estimates, precipitation gauge analyses, and potentially other precipitation estimators at fine time and space scales for the TRMM and GPM eras over the entire globe. The system is run several times for each observation time, first giving a quick estimate and successively providing better estimates as more data arrive. The final step uses monthly gauge data to create research-level products.

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
0.1° - 30 minute|	Gridded, 60°N-60°S, March 2014 to present |	4 months (Prod / final run)	|HDF5 |	FTP (PPS) (/YYYY/MM/DD/3B.HHR.)

----------------------------------------------------------
### radar
----------------------------------------------------------



##### 2B.GPM.DPRGMI.2HCSHv2-1.20150609-S000921-E014155.007256.V04A.HDF5

> 2B-CMB: Combined GMI + DPR single orbit rainfall estimates

> The GPM Combined Radar-Radiometer Algorithm performs two basic functions: first, it provides, in principle, the most accurate, high resolution estimates of surface rainfall rate and precipitation vertical distributions that can be achieved from a spaceborne platform, and it is therefore valuable for applications where information regarding instantaneous storm structure are vital. Second, a global, representative collection of combined algorithm estimates will yield a single common reference dataset that can be used to “cross-calibrate” rain rate estimates from all of the passive microwave radiometers in the GPM constellation. The cross-calibration of radiometer estimates is crucial for developing a consistent, high time-resolution precipitation record for climate science and prediction model validation applications.

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
5 km|	orbital, Past 2 weeks (NRT) |	3 hours (RT); 40 hours (Prod)	|HDF5 |	FTP (PPS) (/YYYY/MM/DD/2B.GPM.DPRGMI)


##### 2A.GPM.Ku.V6-20160118.20150609-S000921-E014155.007256.V04A.HDF5

> The objective of the level 2 DPR algorithms is to generate from the level 1 DPR products radar-only derived meteorological quantities on an instantaneous FOV (field of view) basis. A subset of the results will be used by the level 2 combined radar-radiometer algorithm and the level 3 combined and radar-only products. The general idea behind the algorithms is to determine general characteristics of the precipitation, correct for attenuation and estimate profiles of the precipitation water content, rainfall rate and, when dual-wavelength data are available, information on the particle size distributions in rain and snow. It is particularly important that dual-wavelength data will provide better estimates of rainfall and snowfall rates than the TRMM PR data by using the particle size information and the capability of estimating, even in convective storms, the height at which the precipitation transitions from solid to liquid. 

> 2A-Ku: DPR Ku-only single orbit rainfall estimates

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
5.2km x 125m - 16 orbits per day|	Latitudes 70°N-S, Past 2 Weeks (NRT); March 2014 - present (Prod) |	20 - 120 minutes (NRT); 24 hours (Prod)|HDF5 |	FTP (PPS) (/YYYY/MM/DD/2A.GPM.Ku)

##### 2A.GPM.Ka.V6-20160118.20150609-S014156-E031430.007257.V04A.HDF5

> 2A-Ka: DPR Ka-only single orbit rainfall estimates

Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
5.2km x 125m - 16 orbits per day|	Latitudes 70°N-S, Past 2 Weeks (NRT); March 2014 - present (Prod) |	20 - 120 minutes (NRT); 24 hours (Prod)|HDF5 |	FTP (PPS) (/YYYY/MM/DD/2A.GPM.Ka)



##### 2A.GPM.DPR.V6-20160118.20150609-S000921-E014155.007256.V04A.HDF5

> 2A-DPR: DPR Ka&Ku single orbit rainfall estimates


Resolution | Region - Dates  |	Latency  | 	Format   | 	Source  
-----------|------------------|-----------|-----------|----------
5.2km x 125m - 16 orbits per day|	Latitudes 70°N-S, Past 2 Weeks (NRT); March 2014 - present (Prod) |	20 - 120 minutes (NRT); 24 hours (Prod)|HDF5 |	FTP (PPS) (/YYYY/MM/DD/2A.GPM.DPR)



