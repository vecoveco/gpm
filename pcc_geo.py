
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import wradlib as wrl
import wradlib
from osgeo import osr

'''

    filename = wrl.util.get_wradlib_data_file('geo/Europe_2_02.2017075.aqua.ndvi.250m.tif')
    # pixel_spacing is in output units (lonlat)
    rastercoords, rastervalues = wrl.io.read_raster_data(filename,
                                                         spacing=0.005)
    # specify kwargs for plotting, using terrain colormap and LogNorm
    dem = ax.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                        rastervalues, cmap=plt.cm.terrain, norm=LogNorm(),
                        vmin=1, vmax=3000)
    # make some space on the right for colorbar axis
    div1 = make_axes_locatable(ax)
    cax1 = div1.append_axes("right", size="5%", pad=0.1)
    # add colorbar and title
    # we use LogLocator for colorbar
    cb = plt.gcf().colorbar(dem, cax=cax1,
                           ticks=ticker.LogLocator(subs=range(10)))
    cb.set_label('terrain height [m]')
'''
a = ['20160904','20160917','20160805','20140629','20140921','20141016','20150128','20150427']
xx = range(len(a))

fig = plt.figure(figsize=(20,20))

for i in xx:
    pfad2 = ('/home/velibor/shkgpm/data/'+a[i]+'/dpr/*.HDF5')
    pfad_gpm = glob.glob(pfad2)

    pfad_gpm_g = pfad_gpm[0]

    gpmdpr = h5py.File(pfad_gpm_g, 'r')
    gprof_lat = np.array(gpmdpr['NS']['Latitude'])
    gprof_lon = np.array(gpmdpr['NS']['Longitude'])

    gprof_pp = np.array(gpmdpr['NS']['PRE']['binRealSurface'])#, dtype='float')
    gprof_pp = ((abs((gprof_pp-176.)))*125)

    from pcc import cut_the_swath
    blon, blat, bpp = cut_the_swath(gprof_lon, gprof_lat, gprof_pp, eu=0)

    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    gpm_x, gpm_y = wradlib.georef.reproject(blon, blat, projection_target=proj_stereo , projection_source=proj_wgs)
    grid_xy = np.vstack((gpm_x.ravel(), gpm_y.ravel())).transpose()


    from pcc import plot_borders


    ax = fig.add_subplot(int('33'+str(i+1)), aspect='equal')
    plt.pcolormesh(gpm_x,gpm_y, (bpp*0.125), cmap=plt.cm.terrain, norm=LogNorm(),
                        vmin=1, vmax=3000)
    plot_borders(ax)
    cb = plt.colorbar()
    cb.set_label('binRealSurface [m]')
    plt.grid()
    plt.xlim(-600,400)
    plt.ylim(-4800,-3600)
    plt.title(str(a[i]))

from pcc import plot_dem
ax2 = fig.add_subplot(339, aspect='equal')
#plt.plot(gpm_x[:,0], gpm_y[:,0])
#plt.plot(gpm_x[:,-1], gpm_y[:,-1])
plot_dem(ax2)
plt.title('DEM Height')


plot_borders(ax2)

plt.show()

'''
filename1 = wrl.util.get_wradlib_data_file('geo/radolan_900x900_cr_500.tif')
# pixel_spacing is in output units (lonlat)
rastercoords, rastervalues = wrl.io.read_raster_data(filename1)#,
                                                     #spacing=0.005)



plt.pcolormesh(rastercoords[..., 0], rastercoords[..., 1],
                       rastervalues)
plt.colorbar()
plt.show()
'''
#plt.scatter (h_terra, h_aqua)
#plt.grid()
#plt.xlabel('Terra Hights in m')
#plt.ylabel('Aqua Hights in m')
#plt.show()
