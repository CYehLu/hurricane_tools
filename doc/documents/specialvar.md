# specialvar  

[[source](../.././hurricane_tools//specialvar.py)]  

<span style="color:#a77864">**uv2vrvt_rt**</span>**(u, v, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None)**

    Calculate Vr (radial wind) and Vt (tangential wind) on radius-theta coordinate.
    
    Parameter
    ---------
    u, v : array, shape = (..., ny, nx)
        x and y component wind. Unit is m/s.
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Center longtitude/latitude coordinate.
    radius : 1d array, shape = (nradius,)
        The radius coordinate of data. Unit is km.
    theta : 1d array, shape = (ntheta,). Optional
        The azimuth coordinate of data. Unit is radians.
        It should be "ascent" order (clockwise), and the values should be greater than 0.
    dxdy: Tuple(scalar, scalar). Optional
        Spatial resolution, (dx, dy). Unit is km.
        Default is None, and it would be automatically derived based on `lon` and `lat`.
    intp: str, 'griddata' or 'fortran'. Optional
        See `hurricane_tools.coord_transform.XY2RT` for more information.
        
    Return:
    ------
    vr, vt : array, shape = (..., nradius, ntheta)
        Radial wind and tangential wind on the radius-theta coordinate.
        Unit is m/s.



******
<span style="color:#a77864">**uv2vrvt_xy**</span>**(u, v, lon, lat, clon, clat)**

    Calculate Vr (radial wind) and Vt (tangential wind) on cartesian coordinate.
    
    Parameter:
    ---------
    u, v : n-d array, shape = (..., ny, nx)
        x and y component wind. Unit is m/s
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Center longtitude/latitude coordinate.

    Return:
    ------
    vr, vt : n-d array, shape = (..., ny, nx)
        Radial wind and tangential wind on the cartesian coordinate.
        Unit is m/s.



******
<span style="color:#a77864">**inertial_stability_xy**</span>**(u, v, f, lon, lat, clon, clat, radius=None, thetas=None, dxdy=None, intp=None)**

    Calculate (cyclinic) inertial stability at x-y (longtitude-latitude) coordinate.
    
    Inertial stability is defined as
        I^2 = (f + 2*Vt/r) * (f + 1/r * d(r*Vt)/dr)
    where `f` is coriolis parameter, `Vt` is tangential wind speed, `r` is radius.
    The returned variable is sqrt(I^2).
    
    Parameter:
    ---------
    u, v : n-d array, shape = (..., ny, nx)
        x and y component wind. Unit is m/s.
    f : scalar or array, shape = (ny, nx)
        Coriolis parameter. Unit is s^-1
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Center longtitude/latitude coordinate.
    radius : 1d array, shape = (nradius,). Optional
        The radius coordinate of data. Unit is km.
    theta : 1d array, shape = (ntheta,). Optional
        The azimuth coordinate of data. Unit is radians.
        It should be "ascent" order (clockwise), and the values should be greater than 0.
    dxdy: Tuple(scalar, scalar). Optional
        Spatial resolution, (dx, dy). Unit is km.
        Default is None, and it would be automatically derived based on `lon` and `lat`.
    intp: str, 'griddata' or 'fortran'. Optional
        See `hurricane_tools.coord_transform.XY2RT` for more information.
        
    Return:
    ------
    I : array, shape = (..., ny, nx)
        Inertial stability. Unit is s^-1.



******
<span style="color:#a77864">**inertial_stability_rt**</span>**(vt, f, radius, thetas)**

    Calculate (cyclinic) inertial stability at cylindrical (radius-theta) coordinate.
    
    Inertial stability is defined as
        I^2 = (f + 2*Vt/r) * (f + 1/r * d(r*Vt)/dr)
    where `f` is coriolis parameter, `Vt` is tangential wind speed, `r` is radius.
    
    Parameter:
    ---------
    vt : array, shape = (..., nradius, ntheta)
        Tangential wind speed at cylindrical coordinate.
        Unit is m/s
    f : scalar or 2-d array, shape = (nradius, ntheta)
        Coriolis parameter
        Unit is s^-1
    radius : 1d array, shape = (nradius,)
        Radial coordinate of `vt` and `f`. Unit is km.
    thetas : 1d array, shape = (ntheta,)
        Azimuth coordinate of `vt` and f`. Unit is radian.
        
    Return:
    ------
    I : array, shape = (..., nradius, ntheta)
        Inertial stability. Unit is s^-1.



******