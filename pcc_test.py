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



ZP = '20141007'
year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]

ye = ZP[2:4]

## Read GPM Data
## -------------
pfad2 = ('/home/velibor/shkgpm/data/'+str(year)+str(m)+str(d)+'/dpr/*.HDF5')
pfad_gpm = glob.glob(pfad2)
pfad_gpm_g = pfad_gpm[0]

# GPM Lage und Zeit
gpmdpr = h5py.File(pfad_gpm_g, 'r')
gprof_lat=np.array(gpmdpr['NS']['Latitude'])
gprof_lon=np.array(gpmdpr['NS']['Longitude'])
gpm_time = gpmdpr['NS']['ScanTime']
gpm_zeit = get_time_of_gpm(gprof_lon, gprof_lat, gpm_time)
print 'GPM Datum: ', gpm_zeit


ff = 15
cc = 0.4
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=(12,12))

#GPM Parameter
gpm_para = ['zFactorCorrectedNearSurface','zFactorCorrectedESurface',
            'precipRateNearSurface', 'piaFinal']
gl = len(gpm_para)
for iii in range(gl):
    gprof_pp=np.array(gpmdpr['NS']['SLV'][gpm_para[iii]])
    gprof_pp[gprof_pp==-9999.9]= np.nan


    ## Cut the GPM Swath
    ## ------------------
    blon, blat, gprof_pp_b = cut_the_swath(gprof_lon,gprof_lat,gprof_pp)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)
    gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
    grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()






    ax2 = fig.add_subplot(int('1'+str(gl)+str(iii+1)), aspect='equal')#------------------------------------

    pm2 = plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp_b),
                         cmap=my_cmap, vmin=np.nanmin(gprof_pp_b),
                         vmax=np.nanmax(gprof_pp_b), zorder=2)

    plt.plot(gpm_x[:,0],gpm_y[:,0], color='black',lw=1)
    plt.plot(gpm_x[:,-1],gpm_y[:,-1], color='black',lw=1)
    cb = plt.colorbar(shrink=cc)
    cb.set_label(str(gpm_para[iii]),fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    plt.title('GPM DPR Reflectivity: \n'+ gpm_zeit + ' UTC',fontsize=ff)
    plot_borders(ax2)
    plot_radar(bonnlon, bonnlat, ax2, reproject=True)
    plt.grid(color='r')
    plt.tick_params(
        axis='both',
        which='both',
        bottom='off',
        top='off',
        labelbottom='off',
        right='off',
        left='off',
        labelleft='off')
    plt.xlim(-420,390)
    plt.ylim(-4700, -3700)



plt.tight_layout()
#plt.savefig('/home/velibor/shkgpm/plot/test_dpr_radolan_'+ZP + '.png' )
#plt.close()
plt.show()





