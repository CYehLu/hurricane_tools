import numpy as np
from .distance import latlon2distance
from .transform import uv2vrvt_rt, uv2vrvt_xy
from .fourier import interp_xy_closure


__all__ = [
    'inertial_stability_xy',
    'inertial_stability_rt'
]


def inertial_stability_xy(u, v, f, lon, lat, clon, clat, radius, thetas, dxdy):
    """
    Calculate (cyclinic) inertial stability at x-y (longtitude-latitude) coordinate.
    Inertial stability is defined as
        I^2 = (f + 2*Vt/r) * (f + 1/r * d(r*Vt)/dr)
    where `f` is coriolis parameter, `Vt` is tangential wind speed, `r` is radius.
    
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
    radius : 1d array, shape = (nradius,)
        Radial coordinate (used to calculate the radial gradient)
    thetas : 1d array, shape = (ntheta,)
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
        Inertial stability
    """
    nz, ny, nx = u.shape
    nradius = radius.size
    ntheta = thetas.size

    # calculate vt at radius-theta coordinate, shape = (nz, nradius, ntheta)
    _, vt = uv2vrvt_rt(u, v, lon, lat, clon, clat, radius, thetas, dxdy)
        
    # radial gradient: central difference for interior points, and one-side difference for edge points
    drvt_dr = np.empty((nz, nradius, ntheta))
    rvt = radius[None,:,None] * vt
    drvt_dr[:,1:-1,:] = (rvt[:,2:,:] - rvt[:,:-2,:]) / (radius[2:] - radius[:-2])[None,:,None]
    drvt_dr[:,0,:] = (rvt[:,1,:] - rvt[:,0,:]) / (radius[1] - radius[0])
    drvt_dr[:,-1,:] = (rvt[:,-1,:] - rvt[:,-2,:]) / (radius[-1] - radius[-2])
    
    # interpolate to x-y coordinate
    X = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    Y = latlon2distance(lon, clat, lon, lat)
    X[lon < clon] *= -1
    Y[lat < clat] *= -1
    
    res = np.empty((nz, nx, ny))
    func = interp_xy_closure(radius, thetas, X, Y, center=(0, 0))
    for ilev in range(nz):
        res[ilev,:,:] = func(drvt_dr[ilev,:,:])
    drvt_dr = res
    
    # vt at x-y coordinate
    _, vt = uv2vrvt_xy(u, v, lon, lat, clon, clat)   # (nz, ny, nx)
        
    # finally
    dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    dist_lat = latlon2distance(lon, clat, lon, lat)
    r = 1000 * np.sqrt(X**2 + Y**2)    # km -> m
    r[np.isclose(r, 0)] = np.nan

    I2 = (f + 2*vt/r) * (f + drvt_dr/r)
    I = np.sqrt(np.abs(I2))
    return I


def inertial_stability_rt(vt, f, radius, thetas):
    """
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
    """
    nz, nradius, ntheta = vt.shape
        
    # radial gradient: central difference for interior points, and one-side difference for edge points
    drvt_dr = np.empty((nz, nradius, ntheta))
    rvt = radius[None,:,None] * vt
    drvt_dr[:,1:-1,:] = (rvt[:,2:,:] - rvt[:,:-2,:]) / (radius[2:] - radius[:-2])[None,:,None]
    drvt_dr[:,0,:] = (rvt[:,1,:] - rvt[:,0,:]) / (radius[1] - radius[0])
    drvt_dr[:,-1,:] = (rvt[:,-1,:] - rvt[:,-2,:]) / (radius[-1] - radius[-2])
    
    r = radius[None,:,None].astype(float)
    r[np.isclose(r, 0)] = np.nan
    I2 = (f + 2*vt/r) * (f + drvt_dr/r)
    I = np.sqrt(np.abs(I2))
    return I