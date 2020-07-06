# center.py
ty_center_mslp(lon, lat, slp)
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



******
weighted_ty_center(lon, lat, center_lon, center_lat, var, latlon_range=1)
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



******