import numpy as np
from scipy.optimize import root


__all__ = [
    'latlon2distance',
    'find_dlonlat_by_distance'
]


def latlon2distance(lon1, lat1, lon2, lat2):
    """
    Calculate the distance (km) between two points in longitude/latitude coordinate.
    Units of parameters are degree.
    """
    R = 6373.0    # approximate radius of earth (km)
    
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return R * c


def find_dlonlat_by_distance(clon, clat, distance, xy):
    """
    Find the difference of longitude/latitude which the zonal/meridional
    distance to the (clon, clat) is equal to `distance`.
    
    Parameters
    ----------
    clon, clat: scaler, unit: degree
        The reference coordinate 
    distance: scaler, unit: km
        The distance between target coordinate and (clon, clat).
    xy: str, 'x' or 'y'. 
        If xy == `x`:
            return value = `dlon`, that is the distance between (clon, clat) and
            (clon+dlon, clat) = `distance`.
        If xy == `y`:
            return value = `dlat`, that is the distance between (clon, clat) and
            (clon, clat+dlat) = `distance`.
        
    Return
    ------
    scaler (unit: degree), the difference of longitude or latitude degree.
    """
    if xy == 'x':
        f = lambda ll: latlon2distance(clon, clat, ll, clat) - distance
        return np.abs(root(f, clon).x[0] - clon)
    elif xy == 'y':
        f = lambda ll: latlon2distance(clon, clat, clon, ll) - distance
        return np.abs(root(f, clat).x[0] - clat)