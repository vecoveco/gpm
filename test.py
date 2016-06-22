
import os
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import from_levels_and_colors
import matplotlib.patches as patches
import datetime as dt
import matplotlib.pyplot  as pl
import wradlib
import glob
import h5py
from osgeo import osr

#zdpoly = wradlib.zonalstats.ZonalDataPoly(radar_ll, gprof_gitter[::1], srs=proj_ll, buf=0.)

zd = wradlib.zonalstats.ZonalDataPoly('/user/velibor/SHKGPM/data/test_zonal_poly_cart')
#z = wradlib.zonalstats.ZonalDataPoly.load_vector('/user/velibor/SHKGPM/data/test_zonal_poly_cart/trg.shp')

# GPM
dst_shp = wradlib.util.get_wradlib_data_file('/user/velibor/SHKGPM/data/test_zonal_poly_cart/dst.shp')
d_dataset, d_inLayer = wradlib.io.open_shape(dst_shp)
d_cats, d_keys = wradlib.georef.get_shape_coordinates(d_inLayer)


# Radar
src_shp = wradlib.util.get_wradlib_data_file('/user/velibor/SHKGPM/data/test_zonal_poly_cart/src.shp')
s_dataset, s_inLayer = wradlib.io.open_shape(src_shp)
s_cats, s_keys = wradlib.georef.get_shape_coordinates(s_inLayer)

bbox = s_inLayer.GetExtent()

print(d_cats.shape, d_dataset)
print(s_cats.shape, s_dataset)
print(bbox)


def testplot(cats, catsavg, xy, data,
             levels=[0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 100],
             title=""):
    """Quick test plot layout for this example file
    """
    colors = pl.cm.spectral(np.linspace(0, 1, len(levels)))
    mycmap, mynorm = from_levels_and_colors(levels, colors, extend="max")

    radolevels = [0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 100]
    radocolors = pl.cm.spectral(np.linspace(0, 1, len(radolevels)))
    radocmap, radonorm = from_levels_and_colors(radolevels, radocolors,
                                                extend="max")

    fig = pl.figure(figsize=(10, 16))

    # Average rainfall sum
    ax = fig.add_subplot(211, aspect="equal")
    coll = PatchCollection(cats, array=catsavg, cmap=mycmap, norm=mynorm,
                           edgecolors='white', lw=0.5)
    ax.add_collection(coll)
    ax.autoscale()
    pl.colorbar(coll, ax=ax, shrink=0.5)
    pl.xlabel("GK2 Easting")
    pl.ylabel("GK2 Northing")
    pl.title(title)
    pl.draw()

    # Original radar data
    ax1 = fig.add_subplot(212, aspect="equal")
    pm = pl.pcolormesh(xy[:, :, 0], xy[:, :, 1], np.ma.masked_invalid(data),
                        cmap=radocmap, norm=radonorm)
    coll = PatchCollection(cats, facecolor='None', edgecolor='white', lw=0.5)
    ax1.add_collection(coll)
    cb = pl.colorbar(pm, ax=ax1, shrink=0.5)
    cb.set_label("(mm/h)")
    pl.xlabel("GK2 Easting")
    pl.ylabel("GK2 Northing")
    pl.title("Original radar rain sums")
    pl.draw()
    pl.tight_layout()

testplot(d_cats)