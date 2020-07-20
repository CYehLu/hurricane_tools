# center
<span style="color:#a77864">**tc_center_mslp**</span>**(lon, lat, slp)**

    Find TC center by minimum sea level pressure grid (no weighted).
    
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
<span style="color:#a77864">**weighted_tc_center**</span>**(lon, lat, var, clon=None, clat=None, L=12)**

    Calculate TC center by weighted method.
    
    Parameter:
    ---------
    lon, lat: 
        2-d numpy ndarray. Their shape should be equal.
    var:
        2-d numpy ndarray. Used for determination the center of TC, usually
        is sea level pressure.
        Its shape should equal to lon and lat.
    clon, clat:
        scalar. The first gauess of TC center.
        If None (default), it would use the result of `tc_center_mslp`.
    L:
        int. The half length of weighted box edge length. Default is 12.
        
    Return:
    ------
    tuple, (weighted_center_lon, weighted_center_lat). 



******