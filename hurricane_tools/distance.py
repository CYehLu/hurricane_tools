import numpy as np
from scipy.optimize import root


__all__ = [
    'latlon2distance',
    'find_lonlat_with_distance'
]


def latlon2distance(lon1, lat1, lon2, lat2):
    """calculate the distance (km) of two positions"""
    # approximate radius of earth in km
    R = 6373.0
    
    # convert to radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    return R * c


def find_lonlat_with_distance(clon, clat, distance, xy):
    """
    Find the longtitude/latitude which the horizontal or vertical distance
    to the (clon, clat) is equal to 'distance'.
    
    Parameters:
    ----------
    clon, clat: scaler. The reference coordinate.
    distance: scaler. The distance (km) between target coordinate and (clon, clat).
    xy: str, 'x' or 'y'. 
        If x, it will find the coordinate which distance along latitude line is
        equal to 'distance' (The horizontal distance).
        If y, it will find the coordinate which distance along longtitude line is
        equal to 'distance' (The vertical distance).
        
    Return:
    ------
    scaler, the difference of lon/lat degree
    """
    if xy == 'x':
        f = lambda ll: latlon2distance(clon, clat, ll, clat) - distance
        return np.abs(root(f, clon).x[0] - clon)
    elif xy == 'y':
        f = lambda ll: latlon2distance(clon, clat, clon, ll) - distance
        return np.abs(root(f, clat).x[0] - clat)