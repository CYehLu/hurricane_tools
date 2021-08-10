# center  

[[source](../.././hurricane_tools//center.py)]  

<span style="color:#a77864">**tc_center_mslp**</span>**(lon, lat, slp)**

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



******
<span style="color:#a77864">**weighted_tc_center**</span>**(lon, lat, var, center_fg=None, L=12)**

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



******