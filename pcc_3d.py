"""
============
3D animation
============

A simple example of an animated plot... In 3D!
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")



# draw sphere
u, v = np.mgrid[-180:180:5, -90 : 90:5]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(y, x, z, color="r")


# draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


plt.show()

import math
from wradlib.io import read_generic_netcdf
from wradlib.util import get_wradlib_data_file
import os
import glob

def read_and_overview(filename):
    """Read NetCDF using read_generic_netcdf and print upper level dictionary keys
    """
    test = read_generic_netcdf(filename)
    print("\nPrint keys for file %s" % os.path.basename(filename))
    for key in test.keys():
        print("\t%s" % key)

pfad = ('/automount/ags/velibor/gpmdata/nicht3dComposit/*.nc')

pfad_3d = sorted(glob.glob(pfad))[1]

comp3d = get_wradlib_data_file(pfad_3d)
read_and_overview(comp3d)

from netCDF4 import Dataset
import numpy as np

comp3d = Dataset(pfad_3d, mode='r')



import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Axes3D.plot_trisurf(x,y,z)