# specialvar  

[[source](.././hurricane_tools//specialvar.py)]  

<span style="color:#a77864">**inertial_stability_xy**</span>**(u, v, f, lon, lat, clon, clat, radius=None, thetas=None, dxdy=None)**

    Calculate (cyclinic) inertial stability at x-y (longtitude-latitude) coordinate.
    
    Inertial stability is defined as
        I^2 = (f + 2*Vt/r) * (f + 1/r * d(r*Vt)/dr)
    where `f` is coriolis parameter, `Vt` is tangential wind speed, `r` is radius.
    The returned variable is sqrt(I^2).
    
    Parameter:
    ---------
    u, v : array, shape = (nz, ny, nx)
        Zonal wind and meridional wind at x-y coordinate.
    f : scalar or array, shape = (ny, nx)
        Coriolis parameter
    lon, lat : array, shape = (ny, nx)
        Longtitude and latitude
    clon, clat : scalar
        TC center coordinate
    radius : 1d array, shape = (nradius,). Optional
        Radial coordinate (used to calculate the radial gradient)
    thetas : 1d array, shape = (ntheta,). Optional
        The angles (radians) of each sampled points on the circle.
        See `circular.interp_circle`
        Default is np.arange(*np.deg2rad([0, 360, 1])), the whole circle.
    dxdy: 2-elements tuple, (dx, dy). Optional
        Spatial resolution. 
        Default is None, and it would automatically derive dx and dy besed on `lon`
        and `lat`.
        
    Return:
    ------
    I : array, shape = (nz, ny, nx)
        Inertial stability.



******
<span style="color:#a77864">**inertial_stability_rt**</span>**(vt, f, radius, thetas)**

    Calculate (cyclinic) inertial stability at cylindrical (radius-theta) coordinate.
    
    Inertial stability is defined as
        I^2 = (f + 2*Vt/r) * (f + 1/r * d(r*Vt)/dr)
    where `f` is coriolis parameter, `Vt` is tangential wind speed, `r` is radius.
    
    Parameter:
    ---------
    vt : array, shape = (nz, nradius, ntheta)
        Tangential wind speed at cylindrical coordinate
    f : scalar or array, shape = (nradius, ntheta)
        Coriolis parameter
    radius : 1d array, shape = (nradius,)
        Radial coordinate of `vt` and `f`
    thetas : 1d array, shape = (ntheta,)
        Azimuth coordinate of `vt` and f`
        
    Return:
    ------
    I : array, shape = (nz, nradius, ntheta)
        Inertial stability



******