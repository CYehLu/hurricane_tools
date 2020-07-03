import numpy as np
from .circular import circular_avg_closure


def uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, dxdy=None):
    """u, v: (ny, nx) -> vr, vt: (nradius, ntheta)"""
    
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
    """u, v: (ny, nx) -> vr, vt: (ny, nx)"""
    theta = np.arctan2(lat-clat, lon-clon)   # (ny, nx)
    vr = u * np.cos(theta) + v * np.sin(theta)
    vt = -u * np.sin(theta) + v * np.cos(theta)
    return vr, vt
