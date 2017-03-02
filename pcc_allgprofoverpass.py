"""

This Programm will read the whole GPM GPROF Data from 2014 till 2017
and check der Overpasses with Rain in the Sections of Boxpol and Radolan

"""

import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv

## Read GPM Data
## -------------

pfad = ('/automount/radar-archiv/archiv/GPM/gprof/*.HDF5')
pfad_gprof = glob.glob(pfad)

print 'Es sind ', len(pfad_gprof), ' vorhanden!'

for i in range(10):
    pfad_gprof_g = pfad_gprof[i]

    gpmdprs = h5py.File(pfad_gprof_g, 'r')
    gprof_lat=np.array(gpmdprs['S1']['Latitude'])
    gprof_lon=np.array(gpmdprs['S1']['Longitude'])

    gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])
    gprof_pp[gprof_pp==-9999.9]= np.nan

    print np.nanmax(gprof_pp)

    #with open('/home/velibor/shkgpm/test1.csv', 'wb') as testfile:
    f = open('/home/velibor/shkgpm/test1.csv', 'wb')
    writer = csv.writer(f)
    writer.writerow( ('min', 'max') )
    for i in range(10):
        writer.writerow( (np.nanmin(gprof_pp),(np.nanmin(gprof_pp))) )
    #TODO: Schleife in csv Schreiben lassen
    del(gprof_lat, gprof_lon, gprof_pp)
    #gpm_zeit = get_time_of_gpm(gprof_lon, gprof_lat, gpm_time)
    #print gpm_zeit

f.close()
    #gprof_pp[gprof_pp==-9999.9]= np.nan