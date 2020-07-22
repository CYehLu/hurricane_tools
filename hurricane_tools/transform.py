import numpy as np
from .circular import circular_avg_closure


def uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, dxdy=None):
    """
    Calculate Vr (radial wind) and Vt (tangential wind) on r-theta coordinate (polar coordinate).
    
    Parameter:
    ---------
    u, v : 2d array, shape = (ny, nx)
        x and y component wind.
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Longtitude and latitude of disk center (location of r = 0).
    radius : 1d array, shape = (nradius)
        Radius of circles.
    dxdy : 2 element tuple, (dx, dy). optional
        Spatial resolution. Default is None, and it would find the dx and dy automatically
        based on `lon` and `lat`.
        
    Return:
    ------
    vr, vt : 2d array, shape = (nradius, ntheta)
        Radial wind and tangential wind on the polar coordinate.
    """
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    nradius = radius.size
    
    thetas = np.arange(0, 2*np.pi, 2*np.pi/360)   # (ntheta,)
    
    # interpolate `u` and `v` on the circles
    cirfunc = circular_avg_closure(lon, lat, clon, clat, radius, dxdy, return_interp=True)
    u_interp = cirfunc(u)    # (nradius, ntheta)
    v_interp = cirfunc(v)
    
    vr = u_interp * np.cos(thetas) + v_interp * np.sin(thetas)
    vt = -u_interp * np.sin(thetas) + v_interp * np.cos(thetas)
    return vr, vt


def uv2vrvt_xy(u, v, lon, lat, clon, clat):
    """
    Calculate Vr (radial wind) and Vt (tangential wind) on cartesian coordinate.
    
    Parameter:
    ---------
    u, v : 2d array, shape = (ny, nx)
        x and y component wind.
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude coordinate of `u` and `v`.
    clon, clat : scalar
        Longtitude and latitude of disk center (location of r = 0).

    Return:
    ------
    vr, vt : 2d array, shape = (ny, nx)
        Radial wind and tangential wind on the cartesian coordinate.
    """
    theta = np.arctan2(lat-clat, lon-clon)   # (ny, nx)
    vr = u * np.cos(theta) + v * np.sin(theta)
    vt = -u * np.sin(theta) + v * np.cos(theta)
    return vr, vt
