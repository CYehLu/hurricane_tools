import numpy as np


def ty_center_mslp(lon, lat, slp):
    """
    Find typhoon center by minimum sea level pressure grid (no weighted).
    
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
    
    
def weighted_ty_center(lon, lat, center_lon, center_lat, var, latlon_range=1):
    """
    Giving first guess typhoon center longtitude and latitude, calculate the 
    new center by weighted mothod.
    
    Parameter:
    ---------
    lon, lat: 
        2-d numpy ndarray. Their shape should be equal.
    center_lon, center_lat: 
        float, the first guess typhoon center coordinate.
    var:
        2-d numpy ndarray. Used for determination the center of typhoon, usually
        is pressure.
        Its shape should equal to lon and lat.
    latlon_range: 
        scaler, the degree of lat/lon which will generate a box to calculate the
        weights.
        
    Return:
    ------
    tuple, (new_lon, new_lat). 
    """
    # find the box
    bool_lat = (lat >= center_lat - latlon_range) & (lat <= center_lat + latlon_range)
    bool_lon = (lon >= center_lon - latlon_range) & (lon <= center_lon + latlon_range)
    bool_lonlat = bool_lat & bool_lon
    where = np.where(bool_lonlat)
    left_grid = np.min(where[0])
    right_grid = np.max(where[0])
    bottom_grid = np.min(where[1])
    top_grid = np.max(where[1])
    
    # restrict the lat, lon and var in the box range
    box_lat = lat[left_grid:right_grid,bottom_grid:top_grid]
    box_lon = lon[left_grid:right_grid,bottom_grid:top_grid]
    box_var = var[left_grid:right_grid,bottom_grid:top_grid]
    
    mean_box_var = np.mean(box_var)
    
    # only the positive deviation values can be kept
    weights = np.abs(mean_box_var - box_var) * np.heaviside(mean_box_var - box_var, 0)
    new_center_lon = np.sum(weights * box_lon) / np.sum(weights)
    new_center_lat = np.sum(weights * box_lat) / np.sum(weights)
    
    return (new_center_lon, new_center_lat)