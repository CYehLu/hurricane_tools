import numpy as np
from .circular import interp_circle_closure


__all__ = [
    'uv2vrvt_rt',
    'uv2vrvt_xy'
]


def uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, thetas=None, dxdy=None):
    """
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
    """
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    
    if thetas is None:
        thetas = np.arange(0, 2*np.pi, 2*np.pi/360)   # (ntheta,)
    
    # interpolate `u` and `v` on the circles
    cirfunc = interp_circle_closure(lon, lat, clon, clat, radius, thetas, dxdy)
    
    if u.ndim == 2:
        u_interp = cirfunc(u)    # (nradius, ntheta)
        v_interp = cirfunc(v)
        
    elif u.ndim == 3:
        shape = (u.shape[0], radius.size, thetas.size)
        u_interp = np.empty(shape)
        v_interp = np.empty(shape)
        
        for ilev in range(u.shape[0]):
            u_interp[ilev] = cirfunc(u[ilev,:,:])
            v_interp[ilev] = cirfunc(v[ilev,:,:])
        
    # calc vr/vt at radius/theta coordinate
    vr = u_interp * np.cos(thetas) + v_interp * np.sin(thetas)
    vt = -u_interp * np.sin(thetas) + v_interp * np.cos(thetas)
    return vr, vt


def uv2vrvt_xy(u, v, lon, lat, clon, clat):
    """
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
    """
    theta = np.arctan2(lat-clat, lon-clon)   # (ny, nx)
    vr = u * np.cos(theta) + v * np.sin(theta)
    vt = -u * np.sin(theta) + v * np.cos(theta)
    return vr, vt
