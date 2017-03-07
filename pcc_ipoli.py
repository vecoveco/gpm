#TODO: INTERPOLATION MIT RADIUS!!!
import warnings
from scipy.spatial import cKDTree
import numpy as np
import wradlib.util as util


class MissingSourcesError(Exception):
    """Is raised in case no source coordinates are available for interpolation.
    """
    pass


class MissingTargetsError(Exception):
    """Is raised in case no interpolation targets are available.
    """
    pass

class IpolBase():
    """
    IpolBase(src, trg)

    The base class for interpolation in N dimensions.
    Provides the basic interface for all other classes.

    Parameters
    ----------
    src : ndarray of floats, shape (npoints, ndims)
        Data point coordinates of the source points.
    trg : ndarray of floats, shape (npoints, ndims)
        Data point coordinates of the target points.

    """

    def __init__(self, src, trg):
        src = self._make_coord_arrays(src)
        trg = self._make_coord_arrays(trg)
        self.numsources = len(src)
        self.numtargets = len(trg)

    def __call__(self, vals):
        """
        Evaluate interpolator for values given at the source points.

        Parameters
        ----------
        vals : ndarray of float, shape (numsources, ...)
            Values at the source points which to interpolate

        Returns
        -------
        output : None

        """
        self._check_shape(vals)
        return None

    def _check_shape(self, vals):
        """
        Checks whether the values correspond to the source points

        Parameters
        ----------
        vals : ndarray of float

        """
        assert len(vals) == self.numsources, \
            ('Length of value array %d does not correspond to number '
             'of source points %d' % (len(vals), self.numsources))

    def _make_coord_arrays(self, x):
        """
        Make sure that the coordinates are provided as ndarray
        of shape (numpoints, ndim)

        Parameters
        ----------
        x : ndarray of float with shape (numpoints, ndim)
            OR a sequence of ndarrays of float with len(sequence)==ndim and
            the length of the ndarray corresponding to the number of points

        """
        if type(x) in [list, tuple]:
            x = [item.ravel() for item in x]
            x = np.array(x).transpose()
        elif type(x) == np.ndarray:
            if x.ndim == 1:
                x = x.reshape(-1, 1)
            elif x.ndim == 2:
                pass
            else:
                raise Exception('Cannot deal wih 3-d arrays, yet.')
        return x

class Idw(IpolBase):
    """
    Idw(src, trg, nnearest=4, p=2.)

    Inverse distance weighting interpolation in N dimensions.

    Parameters
    ----------
    src : ndarray of floats, shape (npoints, ndims)
        Data point coordinates of the source points.
    trg : ndarray of floats, shape (npoints, ndims)
        Data point coordinates of the target points.
    nnearest : integer - max. number of neighbours to be considered
    p : float - inverse distance power used in 1/dist**p

    Examples
    --------
    See :ref:`notebooks/interpolation/wradlib_ipol_example.ipynb`.

    Note
    ----
    Uses :class:`scipy:scipy.spatial.cKDTree`

    """

    def __init__(self, src, trg, nnearest=4, p=2.):
        src = self._make_coord_arrays(src)
        trg = self._make_coord_arrays(trg)
        # remember some things
        self.numtargets = len(trg)
        if self.numtargets == 0:
            raise MissingTargetsError
        self.numsources = len(src)
        if self.numsources == 0:
            raise MissingSourcesError
        if nnearest > self.numsources:
            warnings.warn(
                "wradlib.ipol.Idw: <nnearest> is larger than number of "
                "source points and is set to %d corresponding to the "
                "number of source points." % self.numsources,
                UserWarning
            )
            self.nnearest = self.numsources
        else:
            self.nnearest = nnearest
        self.p = p
        # plant a tree
        self.tree = cKDTree(src)
        self.dists, self.ix = self.tree.query(trg, k=self.nnearest)
        # avoid bug, if there is only one neighbor at all
        if self.dists.ndim == 1:
            self.dists = self.dists[:, np.newaxis]
            self.ix = self.ix[:, np.newaxis]
    def __call__(self, vals):

        """
        Evaluate interpolator for values given at the source points.

        Parameters
        ----------
        vals : ndarray of float, shape (numsourcepoints, ...)
            Values at the source points which to interpolate
        maxdist : the maximum distance up to which an interpolated values is
            assigned - if maxdist is exceeded, np.nan will be assigned
            If maxdist==None, values will be assigned everywhere

        Returns
        -------
        output : ndarray of float with shape (numtargetpoints,...)

        """

        # self distances: a list of arrays of distances of the nearest points
        # which are indicated by self.ix
        self._check_shape(vals)
        outshape = list(vals.shape)
        outshape[0] = len(self.dists)
        interpol = (np.repeat(np.nan, util._shape2size(outshape)).
                    reshape(tuple(outshape)).astype('f4'))
        # weights is the container for the weights (a list)
        weights = list(range(len(self.dists)))
        # sources is the container for the source point indices
        src_ix = list(range(len(self.dists)))
        # jinterpol is the jth element of interpol
        jinterpol = 0
        for dist, ix in zip(self.dists, self.ix):
            valid_dists = np.where(np.isfinite(dist))[0]
            dist = dist[valid_dists]
            ix = ix[valid_dists]
            if self.nnearest == 1:
                # defaults to nearest neighbour
                wz = vals[ix]
                w = 1.
            elif dist[0] < 1e-10:
                # if a target point coincides with a source point
                wz = vals[ix[0]]
                w = 1.
            else:
                # weight z values by (1/dist)**p --
                w = 1. / dist ** self.p
                w /= np.sum(w)
                wz = np.dot(w, vals[ix])
            interpol[jinterpol] = wz.ravel()
            weights[jinterpol] = w
            src_ix[jinterpol] = ix
            jinterpol += 1
        return interpol  # if self.qdim > 1  else interpol[0]
