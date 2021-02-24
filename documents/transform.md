# transform  

[[source](.././hurricane_tools//transform.py)]  

<span style="color:#a77864">**uv2vrvt_rt**</span>**(u, v, lon, lat, clon, clat, radius, thetas=None, dxdy=None)**

    Calculate Vr (radial wind) and Vt (tangential wind) on r-theta coordinate (polar coordinate).
    
    Parameter:
    ---------
    u, v : array, shape = (ny, nx) or (nz, ny, nx)
        x and y component wind.
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Longtitude and latitude of disk center (location of r = 0).
    radius : 1d array, shape = (nradius)
        Radius of circles.
    thetas : 1d array, shape = (ntheta)
        The angles (radians) of each sampled points on the circle.
        See `hurricane_tools.circular.interp_circle`
        Default is np.arange(*np.deg2rad([0, 360, 1])), the whole circle.
    dxdy : 2 element tuple, (dx, dy). optional
        Spatial resolution. Default is None, and it would find the dx and dy automatically
        based on `lon` and `lat`.
        
    Return:
    ------
    vr, vt : array, shape = (nradius, ntheta) or (nz, nradius, ntheta)
        Radial wind and tangential wind on the polar coordinate.



******
<span style="color:#a77864">**uv2vrvt_xy**</span>**(u, v, lon, lat, clon, clat)**

    Calculate Vr (radial wind) and Vt (tangential wind) on cartesian coordinate.
    
    Parameter:
    ---------
    u, v : n-d array, shape = (..., ny, nx)
        x and y component wind.
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Longtitude and latitude of disk center (location of r = 0).

    Return:
    ------
    vr, vt : n-d array, shape = (..., ny, nx)
        Radial wind and tangential wind on the cartesian coordinate.



******