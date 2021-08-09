import numpy as np

from .distance import latlon2distance
from .coord_transform import XY2RT, RT2XY
from . import pseudo_coord


__all__ = [
    'uv2vrvt_rt',
    'uv2vrvt_xy',
    'inertial_stability_xy',
    'inertial_stability_rt'
]


def uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None):
    """
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
    """
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    
    if theta is None:
        theta = np.deg2rad(np.arange(360))
    
    xy2rt_obj = XY2RT(lon, lat, clon, clat, radius, theta, dxdy=dxdy, intp=intp)
    
    if u.ndim == 2:
        u_rt = xy2rt_obj(u)    # (nradius, ntheta)
        v_rt = xy2rt_obj(v)

        vr = u_rt * np.cos(theta) + v_rt * np.sin(theta)
        vt = -u_rt * np.sin(theta) + v_rt * np.cos(theta)
        return vr, vt
    
    else:
        shape = u.shape
        ny, nx = shape[-2:]
        u = u.reshape(-1, ny, nx)
        v = v.reshape(-1, ny, nx)
        
        m = u.shape[0]
        u_rt = np.empty((m, radius.size, theta.size))
        v_rt = np.empty((m, radius.size, theta.size))
        vr = np.empty((m, radius.size, theta.size))
        vt = np.empty((m, radius.size, theta.size))
        
        for i in range(m):
            u_rt[i,:,:] = xy2rt_obj(u[i,:,:])
            v_rt[i,:,:] = xy2rt_obj(v[i,:,:])
            
            vr[i,:,:] = u_rt[i,:,:] * np.cos(theta) + v_rt[i,:,:] * np.sin(theta)
            vt[i,:,:] = -u_rt[i,:,:] * np.sin(theta) + v_rt[i,:,:] * np.cos(theta)
            
        vr = vr.reshape(*shape[:-2], radius.size, theta.size)
        vt = vt.reshape(*shape[:-2], radius.size, theta.size)
        return vr, vt
            
    
def uv2vrvt_xy(u, v, lon, lat, clon, clat):
    """
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
    """
    X, Y = pseudo_coord.lonlat2xy(lon, lat, clon, clat)
    theta = np.arctan2(Y, X)
    vr = u * np.cos(theta) + v * np.sin(theta)
    vt = -u * np.sin(theta) + v * np.cos(theta)
    return vr, vt


def inertial_stability_xy(u, v, f, lon, lat, clon, clat, radius=None, thetas=None, dxdy=None, intp=None):
    """
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
    """
    if thetas is None:
        thetas = np.deg2rad(np.arange(360))
        
    if radius is None:
        if dxdy is None:
            dx = latlon2distance(lon[:,1:], lat[:,1:], lon[:,:-1], lat[:,:-1]).mean()
            dy = latlon2distance(lon[1:,:], lat[1:,:], lon[:-1,:], lat[:-1,:]).mean()
            dr = max(dx, dy)
        else:
            dr = max(dxdy)
        
        n = 7
        slc = (slice(n, -n), slice(n, -n))
        dist_x = latlon2distance(clon, lat[slc], lon[slc], lat[slc])
        dist_y = latlon2distance(lon[slc], clat, lon[slc], lat[slc])
        maxdist = min(np.min(dist_x[:,[0,-1]]), np.min(dist_y[[0,-1],:]))
        radius = np.arange(0, maxdist, dr)
        
    shape = u.shape
    ny, nx = shape[-2:]
    nradius = radius.size
    ntheta = thetas.size
    
    # calculate vt at radius-theta coordinate
    # shape = (..., nradius, ntheta) -> (m, nradius, ntheta)
    _, vt = uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, thetas, dxdy)
    vt = vt.reshape(-1, nradius, ntheta)
    m = u.shape[0]
    
    # radial gradient. interior points: central difference, edge points: one-side difference
    rvt = radius[None,:,None] * vt
    drvt_dr = np.empty((m, nradius, ntheta))
    drvt_dr[:,1:-1,:] = (rvt[:,2:,:] - rvt[:,:-2,:]) / (radius[2:] - radius[:-2])[None,:,None]
    drvt_dr[:,0,:] = (rvt[:,1,:] - rvt[:,0,:]) / (radius[1] - radius[0])
    drvt_dr[:,-1,:] = (rvt[:,-1,:] - rvt[:,-2,:]) / (radius[-1] - radius[-2])
    
    # interpolate `drvt_dr` to x-y coordinate. shape = (m, ny, nx)
    tmp = np.empty((m, ny, nx))
    X, Y = pseudo_coord.lonlat2xy(lon, lat, clon, clat)
    rt2xy_obj = RT2XY(lon, lat, clon, clat, radius, thetas, intp=intp)
    for i in range(m):
        tmp[i,:,:] = rt2xy_obj(drvt_dr[i,:,:])
    drvt_dr = tmp
    
    # vt at x-y coordinate
    _, vt = uv2vrvt_xy(u, v, lon, lat, clon, clat)    # (..., ny, nx)
    vt = vt.reshape(m, ny, nx)
    
    # finally
    r = np.sqrt(X**2 + Y**2) * 1000   # km -> m. (ny, nx)
    r[np.isclose(r, 0)] = np.nan
    
    I2 = (f + 2*vt/r) * (f + drvt_dr/r)
    I = np.sqrt(np.abs(I2))
    I = I.reshape(shape)
    return I


def inertial_stability_rt(vt, f, radius, thetas):
    """
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
    """
    shape = vt.shape
    nradius, ntheta = shape[-2:]
    
    vt = vt.reshape(-1, nradius, ntheta)
    m = vt.shape[0]
    
    radius = radius * 1000   # km -> m
        
    # radial gradient. interior points: central difference, edge points: one-side difference
    rvt = radius[None,:,None] * vt
    drvt_dr = np.empty((m, nradius, ntheta))
    drvt_dr[:,1:-1,:] = (rvt[:,2:,:] - rvt[:,:-2,:]) / (radius[2:] - radius[:-2])[None,:,None]
    drvt_dr[:,0,:] = (rvt[:,1,:] - rvt[:,0,:]) / (radius[1] - radius[0])
    drvt_dr[:,-1,:] = (rvt[:,-1,:] - rvt[:,-2,:]) / (radius[-1] - radius[-2])
    
    r = radius[None,:,None].astype(float)
    r[np.isclose(r, 0)] = np.nan
    I2 = (f + 2*vt/r) * (f + drvt_dr/r)
    I = np.sqrt(np.abs(I2))
    I = I.reshape(shape)
    return I