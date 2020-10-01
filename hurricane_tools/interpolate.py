import numpy as np
from scipy.spatial import qhull


class FastGriddata:
    """
    Faster version of scipy.interpolate.griddata (in repeatly interpolating 
    case).
    
    When performs interpolation repeatly, and parameters remain the same
    except the interpolated values, this can reduce some computational 
    time by reuse the "triangle grid" during the process of griddata.
    
    Codes are modified from : https://reurl.cc/xZ50Yb
    The first and the second answers.
    
    Example:
    -------
    >>> import time
    >>> from scipy.interpolate import griddata
    >>> from hurricane_tools.interpolate import FastGriddata
    >>>
    >>> # construct data points
    ... X_, Y_ = np.meshgrid(np.linspace(-10, 10), np.linspace(-10, 10))
    ... xy = np.vstack((X_.ravel(), Y_.ravel())).T   # shape=(100*100, 2)
    ... uv = np.random.rand(1000, 2) * 20 - 10
    >>>
    >>> # construct 3 different interpolating values
    ... vals1 = np.ravel(np.cos(X_) * np.cos(Y_))
    ... vals2 = np.ravel(np.sin(X_) + 2 * np.cos(Y_))
    ... vals3 = np.ravel(np.cos(X_**2) * np.cos(Y_))
    >>>
    >>> # using scipy.interpolate.gridata -- it is slower
    ... time_start = time.time()
    ... res1_g = griddata(xy, vals1, uv)
    ... res2_g = griddata(xy, vals2, uv)
    ... res3_g = griddata(xy, vals3, uv)
    ... time_end = time.time()
    ... print('by using griddata :', time_end - time_start)   # 0.1519 sec
    by using griddata : 0.15190553665161133
    >>>
    >>> # using FasterGriddata -- faster than apply griddata 3 times
    ... time_start = time.time()
    ... interp_obj = FastGriddata(xy, uv)
    ... res1_f = interp_obj.interpolate(vals1)
    ... res2_f = interp_obj.interpolate(vals2)
    ... res3_f = interp_obj.interpolate(vals3)
    ... time_end = time.time()
    ... print('by using FastGriddata :', time_end - time_start)   # 0.0543 sec
    by using FastGriddata : 0.05439400672912598
    >>> 
    >>> # check results
    ... np.allclose(res1_g, res1_f)
    True
    >>> np.allclose(res2_g, res2_f)
    True
    >>> np.allclose(res3_g, res3_f)
    True
    """
    
    def __init__(self, xy, uv, d=2, fill_value=np.nan):
        """
        Initialization for interpolation.
        
        Parameter
        ---------
        xy : 2-d array, shape = (n, 2)
            Coordinate of `values` in `FastGriddata.interpolate`
        uv : 2-d array, shape = (m, 2)
            Coordinate of interpolation result.
        d : int
            Default is 2.
        fill_value : scalar
            Value used to fill in for the extrapolation points.
            Default is `nan`.
        """
        tri = qhull.Delaunay(xy)
        simplex = tri.find_simplex(uv)       
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uv - temp[:,d]
        bary = np.einsum('njk,nk->nj', temp[:,:d,:], delta)
        
        # calculate interpolation weights, and dealing with extrapolation points
        weights = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))
        weights[simplex == -1, :] = fill_value
        
        self.vertices = vertices
        self.weights = weights
    
    def interpolate(self, values):
        """
        Interpolate `values` from `xy` coordinate to `uv` coordinate.
        
        Parameter
        ---------
        values : 2-d array, shape = (n, 2)
        
        Return
        ------
        2-d array, shape = (m, 2)
        """
        return np.einsum('nj,nj->n', np.take(values, self.vertices), self.weights)