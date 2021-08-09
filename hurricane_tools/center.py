import numpy as np


__all__ = [
    'tc_center_mslp',
    'weighted_tc_center'
]


def tc_center_mslp(lon, lat, slp):
    """
    Find TC center by finding the grid of minimum sea level pressure.
    
    Parameter
    ---------
    lon, lat : 2d numpy array, shape = (ny, nx)
        longtitude / latitude
    slp : 2d numpy array, shape = (ny, nx)
        sea level pressure
        
    Return
    ------
    tuple, (center_lon, center_lat)
    """
    i, j = np.unravel_index(np.argmin(slp), slp.shape)
    clon, clat = lon[i,j], lat[i,j]
    return clon, clat
    

def weighted_tc_center(lon, lat, var, center_fg=None, L=12):
    """
    Calculate TC center by weighted method.
    
    Parameter:
    ---------
    lon, lat: 
        2-d numpy ndarray. Their shape should be equal.
    var:
        2-d numpy ndarray. Used for determination the center of TC, usually
        is sea level pressure.
        Its shape should equal to lon and lat.
    center_fg:
        Tuple(scalar, scalar). Optional
        The first guess of TC center (lon, lat).
        If None, it would use the result of `tc_center_mslp`.
    L:
        int. The half length of weighted box edge length. Default is 12.
        
    Return:
    ------
    Tuple(scalar, scalar)
        weighted TC center (lon, lat). 
    """
    if center_fg is None:
        clon, clat = tc_center_mslp(lon, lat, var)
    else:
        clon, clat = center_fg
    
    # find the nearest grid point to clon and clat
    diff_lon = np.abs(lon - clon)
    diff_lat = np.abs(lat - clat)
    idx = np.unravel_index(np.argmin(diff_lon + diff_lat), lon.shape)
    
    # box area
    lon_b = lon[idx[0]-L:idx[0]+L, idx[1]-L:idx[1]+L]
    lat_b = lat[idx[0]-L:idx[0]+L, idx[1]-L:idx[1]+L]
    var_b = var[idx[0]-L:idx[0]+L, idx[1]-L:idx[1]+L]
    
    # find weighted center
    weight = np.zeros_like(var_b)
    mean_var_b = np.mean(var_b)
    wloc = var_b < mean_var_b
    weight[wloc] = (mean_var_b - var_b)[wloc]
    sumw = np.sum(weight)
    new_clon = np.sum(weight * lon_b) / sumw
    new_clat = np.sum(weight * lat_b) / sumw
    return new_clon, new_clat