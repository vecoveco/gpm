"""

This Programm will read the whole GPM GPROF Data from 2014 till 2017
and check der Overpasses with Rain in the Sections of Boxpol and Radolan

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
from pcc import cut_the_swath2

## Landgrenzenfunktion
## -------------------
from pcc import boxpol_pos
bonn_pos = boxpol_pos()
bx, by = bonn_pos['gkx_ppi'], bonn_pos['gky_ppi']
bonnlat, bonnlon = bonn_pos['lat_ppi'], bonn_pos['lon_ppi']
from pcc import plot_borders
from pcc import plot_radar
#from pcc import get_miub_cmap
#my_cmap = get_miub_cmap()

from pcc import get_my_cmap
my_cmap = get_my_cmap()
from pcc import cut_the_swath

GGG = []
RRR = []

maxi, mini = [], []
import csv
## Read GPM Data
## -------------

pfad = ('/automount/radar-archiv/archiv/GPM/gprof/*.HDF5')
pfad_gprof = sorted(glob.glob(pfad))


bonn_lat1 = 46.952580411190304
bonn_lat2 = 54.896591448461479
bonn_lon1 = 2.0735617005681379
bonn_lon2 = 15.704155593113517
r1 = np.array([bonn_lon1, bonn_lon2, bonn_lon2, bonn_lon1, bonn_lon1])
r2 = np.array([bonn_lat1, bonn_lat1, bonn_lat2, bonn_lat2, bonn_lat1])

print 'Es sind ', len(pfad_gprof), ' vorhanden!'



f = open('output3bonndpr.csv','w')

for i in range(len(pfad_gprof)):
    pfad_gprof_g = pfad_gprof[i]
    print i, ' von ', len(pfad_gprof)
    #print pfad_gprof_g
    try:
        gpmdprs = h5py.File(pfad_gprof_g, 'r')
        #[:,86:135]
        gprof_lat=np.array(gpmdprs['S1']['Latitude'])#[:,86:135]
        gprof_lon=np.array(gpmdprs['S1']['Longitude'])#[:,86:135]
        gprof_pp=np.array(gpmdprs['S1']['surfacePrecipitation'])#[:,86:135]
        gpm_time = gpmdprs['S1']['ScanTime']#[:,86:135]
        #gpm_zeit = get_time_of_gpm(gprof_lon[0,0], gprof_lat[0,0], gpm_time[0,0])
        #print gpm_zeit
        gprof_pp[gprof_pp<0]= np.nan

        gprof_lat, gprof_lon, gprof_pp = gprof_lat[:,86:135], gprof_lon[:,86:135], gprof_pp[:,86:135]


        from pcc import cut_the_swath
        blat, blon, bpp = cut_the_swath(gprof_lat, gprof_lon, gprof_pp, eu=2)


        #with open("output1.csv",'wb') as f:
        try:
            if np.nansum(bpp) >= 0:
                writer = csv.writer(f, dialect='excel')
                writer.writerow([pfad_gprof_g[-41::],np.nanmax(bpp), np.nanmin(bpp), np.nanmean(bpp)])

                #plt.subplot(1,2,1)
                #plt.pcolormesh(blat,blon,np.ma.masked_invalid(bpp))
                #plt.plot(r1,r2)
                #plt.xlim(-10,26)
                #plt.ylim(35,65)
                #plt.subplot(1,2,2)
                #plt.pcolormesh(blat[:,86:135], blon[:,86:135], bpp[:,86:135])
                #plt.plot(r1,r2)
                #plt.xlim(-10,26)
                #plt.ylim(35,65)
                #plt.show()
        except:
            pass

        del(pfad_gprof_g, gpmdprs, gprof_lon, gprof_lat, gprof_pp, gpm_time,blat, blon, bpp)
    except:
        pass

f.close()
    #



from pcc import melde_dich
melde_dich('Das Program RADOLAN mit DPR ist fertig!')
'''
G_all = np.concatenate(GGG,axis=0)
R_all = np.concatenate(RRR,axis=0)
from pcc import plot_scatter


fig = plt.figure(figsize=(12,12))
ax11 = fig.add_subplot(111, aspect='equal')
plot_scatter(G_all, R_all)
import matplotlib as mpl
mean = [ np.nanmean(G_all),np.nanmean(R_all)]
width = np.nanstd(G_all)
height = np.nanstd(R_all)
angle = 0
ell = mpl.patches.Ellipse(xy=mean, width=width, height=height,
                          angle=180+angle, color='blue', alpha=0.8,
                          fill=False, ls='--', label='Std')
ax11.add_patch(ell)
plt.xlabel(('GPM DPR (dBZ)'))
plt.ylabel(('RADOLAN (dBZ)'))
plt.grid()
plt.savefig('/home/velibor/shkgpm/plot/allegprof/all_gpm_gprof_radolan_'+ r_pro  + '.png' )
plt.close()
'''



# Als Datei abspeichern
#import pandas as pd
#df = pd.DataFrame({'max':maxi,'min':mini})
#print df.shape
#df.to_csv("/home/velibor/shkgpm/test1.csv", sep=',')