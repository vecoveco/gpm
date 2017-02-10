
import h5py
import numpy as np
import matplotlib.pyplot as plt
import wradlib
import glob
from scipy import stats
import wradlib as wrl
from osgeo import osr


#ZP = '20141007023500'
ZP = '20161024232500'
year, m, d, ht, mt, st = ZP[0:4], ZP[4:6], ZP[6:8], ZP[8:10], ZP[10:12], ZP[12:14]
ye = ZP[2:4]

## Read RADOLAN GK Koordinaten
## ----------------------------

pfad = ('/automount/radar/dwd/rx/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/raa01-rx_10000-'+str(ye)+str(m)+
        str(d)+str(ht)+str(mt)+'-dwd---bin.gz')

pfad_radolan = pfad[:-3]

rw_filename = wradlib.util.get_wradlib_data_file(pfad)
rwdata, rwattrs = wradlib.io.read_RADOLAN_composite(rw_filename)


radolan_grid_xy = wradlib.georef.get_radolan_grid(900,900)
x = radolan_grid_xy[:,:,0]
y = radolan_grid_xy[:,:,1]


rx = x[np.where(rwdata != -9999)]
ry = y[np.where(rwdata != -9999)]
rd = rwdata[np.where(rwdata != -9999)]

nans = rwdata == -9999

rwdata = np.ma.masked_equal(rwdata, -9999) / 2 - 32.5


gk3 = wradlib.georef.epsg_to_osr(31467)
proj_stereo = wrl.georef.create_osr("dwd-radolan")
proj_wgs = osr.SpatialReference()
proj_wgs.ImportFromEPSG(4326)


# Radolan Shape file erstellen
rrrr = np.vstack((rx.ravel(), ry.ravel())).transpose()
rado = wradlib.zonalstats.DataSource(rrrr, srs=gk3)
rado.dump_vector('rado.shp')

xn = x.copy()
yn = y.copy()


xn[nans] = np.nan
yn[nans] = np.nan

firstx = np.nanargmax(xn > -600, axis=1)
lastx = np.nanargmax(xn[:,::-1] < 500, axis=1)


firsty = np.nanargmax(yn > -5000., axis=0)
lasty = np.nanargmax(yn[::-1,:] < -3600., axis=0)


plt.subplot(3,2,1)
plt.plot(firstx)
plt.title('first x')
plt.subplot(3,2,2)
plt.plot(firsty)
plt.title('first y')
plt.subplot(3,2,3)
plt.plot(lastx)
plt.title('last x')
plt.subplot(3,2,4)
plt.plot(lasty)
plt.title('last y')
plt.subplot(3,2,5)
plt.imshow(xn)
plt.show()


#Radolanrand
def cut_the_edge(gpm_x, gpm_y, rrr):
    rand_y_unten = -4658.6447242655722
    rand_y_oben = -3759.6447242655722
    rand_x_rechts = 375.5378330781441


    gy = gpm_y.copy()
    gx = gpm_x.copy()
    gg = rrr.copy()

    gg[np.where(gy < rand_y_unten)] = np.nan
    gg[np.where(gy > rand_y_oben)] = np.nan
    gg[np.where(gx > rand_x_rechts)] = np.nan


'''
# https://github.com/dwyerk/boundaries/blob/master/concave_hulls.ipynb
# alpha shape funktion
import fiona
import shapely.geometry as geometry
from descartes import PolygonPatch

def plot_polygon(polygon):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    margin = .3
    x_min, y_min, x_max, y_max = polygon.bounds
    ax.set_xlim([x_min-margin, x_max+margin])
    ax.set_ylim([y_min-margin, y_max+margin])
    patch = PolygonPatch(polygon, fc='#999999',
                         ec='#000000', fill=True,
                         zorder=-1)
    ax.add_patch(patch)
    return fig

from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Delaunay
import numpy as np
import math



def alpha_shape(points, alpha):

    if len(points) < 4:
        # When you have a triangle, there is no sense in computing an alpha
        # shape.
        return geometry.MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):

        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])

    coords = np.array([point.coords[0] for point in points])

    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]

        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)

        # Semiperimeter of triangle
        s = (a + b + c)/2.0

        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)

        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points


import fiona
import shapely.geometry as geometry
input_shapefile = '/automount/user/velibor/SHKGPM/prog/rado.shp'
shapefile = fiona.open(input_shapefile)
points = [geometry.shape(point['geometry'])
          for point in shapefile]

concave_hull, edge_points = alpha_shape(points,alpha=0.7)
plot_polygon(concave_hull)
plt.plot(rx,ry,'o', color='#f16824')
plt.show()



rado_concave_hull = wradlib.zonalstats.DataSource(concave_hull, srs=gk3)
rado.dump_vector('datapro/rado_concave_hull.shp')

rado_edge_points = wradlib.zonalstats.DataSource(edge_points, srs=gk3)
rado.dump_vector('datapro/rado_edge_points.shp')

zd = wradlib.zonalstats.ZonalDataPoly('/user/velibor/SHKGPM/prog/datapro/rado_concave_hull.shp')
rado_hull= wradlib.zonalstats.ZonalDataPoly.load_vector('/user/velibor/SHKGPM/prog/datapro/rado_concave_hull.shp')
'''