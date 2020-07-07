# circular.py
<span style="color:#a77864">**circular_avg**</span>**(lon, lat, values, clon, clat, radius, dxdy=None, **kwargs)**

    Calculate circular mean.
    
    Parameters:
    ----------
    lon, lat: 2-d array, shape=(ny, nx)
        The lon/lat coordinates of `values`.
    values: 2-d array, shape=(ny, nx)
        Data values
    clon, clat: scaler
        The center coordinates of circular mean.
    radius: scaler or 1-d array-like
        The radius of circles. Unit is km.
    dxdy: 2-elements tuple, (dx, dy). Optional
        Spatial resolution. 
        Default is None, and it would automatically derive dx and dy besed on `lon`
        and `lat`.
        
    **kwargs: {theta, return_interp, method, fill_value, rescale}
        theta: 3-elements tuple. (theta_start, theta_end, dtheta)
            The start angle on the circle (theta_start), the ending angle on the
            circle (theta_end), and the angle interval of the circle samples (dtheta).
            Default is (0, 2*pi, 2*pi/360), the whole circle.
        return_interp: bool. Optional
            If `return_interp` = True, it would return the whole interpolating result
            with shape = (n_radius, n_theta), where `n_radius` is the number of radius
            samples and `n_theta` is the number of theta samples.
            And the circular mean is just the theta average of return result.
            e.g >>> interp = circular_avg(..., return_interp=True)   # shape=(n_radius, n_theta)
                >>> cir_avg = interp.mean(axis=1)   # shape=(n_radius,)
            Default is False.
        method, fill_value, rescale:
            The parameters used in scipy.interpolate.griddata.
            Check document for griddata:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
            
    Returns:
    -------
    circular average result, 1-d array with size=len(radius).
        
    NOTE:
    ----
    If this function is used repeatly, and all parameters remain the same except `values`, it is 
    worth to use `circular_avg_closure` instead. `circular_avg_closure` returns a closure function
    which only use `values` as its parameter.



******
<span style="color:#a77864">**circular_avg_closure**</span>**(lon, lat, clon, clat, radius, dxdy=None, **kwargs)**

    Calculate circular mean.
    
    This function is very similar to `circular_avg`, while this function returns
    a closure function which use `values` as its parameters.
    This is suitable for the situations which needed to calculate circular mean
    repeatly, and only `vaules` are different, all other parameters remain the same.
    
    Parameters:
    ----------
    See `circular_avg`
    
    Returs:
    ------
    closure function



******
<span style="color:#a77864">**rmw**</span>**(lon, lat, ws, clon, clat, maxdist=550, box=True, **kwargs)**

    Find TC RMW
    
    Paramters
    ---------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude / latitude
    ws : 2d array, shape = (ny, nx)
        Wind speed
    clon, clat : scalar
        TC center longtitude / latitude
    maxdist : scalar
        The maximum search distance (km). default is 550
    box : bool
        Create a box area that only calculate interpolation in this
        box. Default is True.
    **kwargs : 
        Keyword arguments for `circular_avg`
        
    Return
    ------
    dist_coord : 1d array, shape = (n,)
        The coordinate of `axissym_ws`.
        The distance between i'th `axissym_ws` profile location
        and TC center is dist_coord[i].
    axissym_ws : 1d array, shape = (n,)
        Axis-symmetric wind speed profile
    rmw : scalar
        Radius of maximum wind speed



******
<span style="color:#a77864">**axisymmetricity**</span>**(lon, lat, var, radius, clon, clat, dxdy=None, integ='trapz', **kwargs)**

    Calculate axisymmetricity based on Miyamoto and Takemi (2013).
    
    Parameter:
    ---------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude / latitude
    var : 2d array, shape = (ny, nx)
        The variable used in the calculation of axissymmetricity, e.g: wind speed.
    radius : 1d array, shape = (nr,)
        The radius of circles. Unit is km.
    clon, clat: scaler
        The center longtitude / latitude of TC.
    dxdy : 2 element tuple, (dx, dy).
        Spatial resolution. Default is None, and it would find the dx and dy automatically
        based on `lon` and `lat`.
    integ : str, {'trapz', 'simps'}.
        Numerical integration method. 'trapz' is trapezoidal method, and `simps` is Simpson’s
        method. 
        See scipy document: https://reurl.cc/X6KpYD
    **kwargs : method, fill_value, rescale
        The parameters used in scipy.interpolate.griddata.
        Check document for griddata:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
        
    Return:
    ------
    1d array with shape = (nr,). The axisymmetricity at the given radius.
    
    Reference
    ---------
    [1] Yoshiaki Miyamoto and Tetsuya Takemi: "A Transition Mechanism for the Spontaneous 
        Axisymmetric Intensification of Tropical Cyclones"
        J. Atmos. Sci, 70, 112-129
        https://doi.org/10.1175/JAS-D-11-0285.1



******