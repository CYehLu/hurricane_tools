import numpy as np
from .distance import latlon2distance, find_lonlat_with_distance


__all__ = [
    'lonlat2xy',
    'xy2lonlat'
]


def lonlat2xy(lon, lat, clon, clat):
    """
    Convert to distance-based coordinate.
    
    Parameter:
    ---------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude
    clon, clat : scalar
        Center longtitude-latitude coordinate
        
    Return:
    ------
    X, Y : 2d array, shape = (ny, nx)
        `X` is the zonal distance between each lon-lat coordinate points and
        the center coordinate, and `Y` is the meridional distance.
    """
    X = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    Y = latlon2distance(lon, clat, lon, lat)
    X[lon < clon] *= -1
    Y[lat < clat] *= -1
    return X, Y


def xy2lonlat(X, Y, clon=None, clat=None, equal_dist=True):
    """
    Convert to longtitude-latitude coordinate.
    
    Parameter:
    ---------
    X, Y : 2d array, shape = (ny, nx)
        `X` is the zonal distance between each grid points and the center
        coordinate (`clon`, `clat`), and `Y` is the meridional distance.
    clon, clat : scalar, optional
        Center longtitude-latitude coordinate.
        Default is clon=120, clat=20
    equal_dist : bool, optional
        Determine if the zonal/meridional distance between each zonal/meridional
        grids is nearly equal.
        This can only be True now.
        
    Return:
    ------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and Latitude.
    """
    clon = 120 if clon is None else clon
    clat = 20 if clat is None else clat
    
    ny, nx = X.shape
    
    if equal_dist:
        dx = np.diff(X, axis=1).mean(axis=1)   # (ny,)
        dy = np.diff(Y, axis=0).mean(axis=0)   # (nx,)
    else:
        raise ValueError("not finished...")
        
    # calculate `dlon` (shape=(ny,)) and `dlat` (shape=(nx,))
    if np.allclose(dx, dx[0]):
        dlon = np.ones((ny,)) * find_lonlat_with_distance(clon, clat, dx[0], 'x')
    else:
        dlon = np.array([find_lonlat_with_distance(clon, clat, d, 'x') for d in dx])
        
    if np.allclose(dy, dy[0]):
        dlat = np.ones((nx,)) * find_lonlat_with_distance(clon, clat, dy[0], 'y')
    else:
        dlat = np.array([find_lonlat_with_distance(clon, clat, d, 'y') for d in dy])
        
    isevenx = (nx+1) % 2
    iseveny = (ny+1) % 2
    lon = np.linspace(clon-dlon*(nx//2), clon+dlon*(nx//2-isevenx), nx).T   # (ny, nx)
    lat = np.linspace(clat-dlat*(ny//2), clat+dlat*(ny//2-iseveny), ny)     # (ny, nx)
    return lon, lat