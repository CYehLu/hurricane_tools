import numpy as np
from scipy.interpolate import griddata
from scipy.integrate import trapz, simps

from .distance import latlon2distance
from .interpolate import FastGriddata


def _interp_circle_xy(X, Y, values, cx, cy, radius, theta, dxdy):
    """interpolating data on the circle for cartesain coordinate"""
    raise ValueError("This function is not finished yet")


def _interp_circle_lonlat(lon, lat, values, clon, clat, radius, theta, dxdy, **kwargs):
    """interpolating data on the circle for lon/lat coordinate"""
    dx, dy = dxdy
    nums_r = radius.size
    nums_t = theta.size
    
    # the meridional/zonal distance between center and every grid points
    dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    dist_lat = latlon2distance(lon, clat, lon, lat)
    dist_full = np.sqrt(dist_lon**2 + dist_lat**2)
    
    # give signs to dist_lon/lat (west of `clon` or south of `clat` would be negative)
    dist_lon[lon < clon] *= -1
    dist_lat[lat < clat] *= -1
    
    # set a box area to reduce the amount of computation
    max_r = radius.max()
    L_lon = int(max_r // dx) + 6
    L_lat = int(max_r // dy) + 6
    cix = np.unravel_index(dist_full.argmin(), dist_full.shape)  # (nearly) center index
    dist_lon_b = dist_lon[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
    dist_lat_b = dist_lat[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
    values_b = values[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
    
    # reshape into interpolation form
    dist_lonlat_b = np.vstack((dist_lon_b.ravel(), dist_lat_b.ravel())).T
    val_b = values_b.ravel()
    
    # construct circular samples points and interpolate
    circle_pts = np.zeros((nums_t*nums_r, 2))
    circle_pts[:,0] = radius.repeat(nums_t) * np.cos(np.tile(theta, nums_r))
    circle_pts[:,1] = radius.repeat(nums_t) * np.sin(np.tile(theta, nums_r))
    interp = griddata(dist_lonlat_b, val_b, circle_pts, **kwargs)
    
    return interp.reshape(nums_r, nums_t)
    

def interp_circle(X, Y, values, cx, cy, radius, theta=None, dxdy=None, coord='lonlat', **kwargs):
    """
    Interpolating data on the circles.
    
    Parameters:
    ----------
    X, Y: 2-d array, shape = (ny, nx)
        The coordinates of the `values`.
    values: 2-d array, shape = (ny, nx)
        Values on `X` and `Y` coordinate.
    cx, cy: scalar
        Center coordinate of the circle.
    radius: scalar or 1-d array-like, shape = (nradius,)
        The radius of circles.
    theta: 1-d array, shape = (ntheta,). Optional
        The angles (radians) of each sampled points on the circle.
        Default is np.arange(*np.deg2rad([0, 360, 1])), the whole circle.
        Examples:
        >>> theta = np.arange(*np.deg2rad([0, 360, 1]))   # whole circle
        >>> theta = np.arange(*np.deg2rad([0, 180, 1]))   # upper half circle
        >>> theta = np.arange(*np.deg2rad([90, 270, 10]))   # left half circle, coarser samples
    dxdy: 2-elements tuple, (dx, dy). Optional
        Spatial resolution. 
        Default is None, and it would automatically derived based on `X` and `Y`.
    coord: str, 'lonlat' or 'xy'. Optional
        The coordinate system of `X` and `Y`.
        If coord = `xy`, then `X` and `Y` are cartesain coordinate.
        If coord = 'lonlat', then `X` and `Y` are longtitude and latitude.
        Default is `lonlat`.
    **kwargs: {method, fill_value, rescale}
        The parameters used in scipy.interpolate.griddata.
        Check document for griddata:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
        
    Return:
    ------
    Interpolating result, shape = (nradius, ntheta)
    """
    # convert `radius` to iterable
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    
    if theta is None:
        theta = np.arange(*np.deg2rad([0, 360, 1]))
    
    if coord == 'xy':
        if dxdy is None:
            dx = X[0,1] - X[0,0]
            dy = Y[1,0] - Y[0,0]
            dxdy = (dx, dy)
        return _interp_circle_xy(X, Y, values, cx, cy, radius, theta, dxdy, **kwargs)
    
    elif coord == 'lonlat':
        if dxdy is None:
            dx = latlon2distance(X[0,0], Y[0,0], X[0,1], Y[0,0])
            dy = latlon2distance(X[0,0], Y[0,0], X[0,0], Y[1,0])
            dxdy = (dx, dy)
        return _interp_circle_lonlat(X, Y, values, cx, cy, radius, theta, dxdy, **kwargs)
    
    else:
        raise ValueError(f"Unavailable coord: {coord}. It shold be 'lonlat' or 'xy'.")


def circular_avg(lon, lat, values, clon, clat, radius, dxdy=None, **kwargs):
    """
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
        
    **kwargs: {theta, method, fill_value, rescale}
        theta: 1-d array
            The angles (radians) of each sampled points on the circle.
            See `interp_circle`.
            Default is np.arange(*np.deg2rad([0, 360, 1])), the whole circle.
        method, fill_value, rescale:
            The parameters used in scipy.interpolate.griddata.
            Check document for griddata:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
            
    Returns:
    -------
    The circular average result on each radius, shape = (len(radius),)
        
    NOTE:
    ----
    If this function is used repeatly, and all parameters remain the same except `values`, it is 
    worth to use `circular_avg_closure` instead. `circular_avg_closure` returns a closure function
    which only use `values` as its parameter.
    """
    theta = kwargs.pop('theta', None)
    res = interp_circle(lon, lat, values, clon, clat, radius, theta, dxdy, 'lonlat', **kwargs)
    return res.mean(axis=1)


def circular_avg_closure(lon, lat, clon, clat, radius, dxdy=None, **kwargs):
    """
    Return a closure function to calculate circular mean.
    
    This function is very similar to `circular_avg`, while this function returns
    a closure function which use `values` as its parameters.
    This is suitable for the situations which needed to calculate circular mean
    repeatly, and only `vaules` are different, all other parameters remain the same.
    
    Parameters:
    ----------
    See `circular_avg`
    
    Returs:
    ------
    A closure function, which its parameter is `values` and return the interpolateing
    result (see `circular_avg` for `values`).
    
    Example:
    -------
    >>> # using `circular_avg`
    >>> res1 = circular_avg(lon, lat, val1, clon, clat, radius)
    >>> res2 = circular_avg(lon, lat, val2, clon, clat, radius)
    >>> 
    >>> # using `circular_avg_closure`
    >>> cavg_func = circular_avg_closure(lon, lat, clon, clat, radius)
    >>> res1_closure = cavg_func(val1)
    >>> res2_closure = cavg_func(val2)
    >>> 
    >>> np.allclose(res1, res1_closure)
    True
    >>> np.allclose(res2, res2_closure)
    True
    """
    # get `theta`
    theta = kwargs.pop('theta', None)
    if theta is None:
        theta = np.arange(*np.deg2rad([0, 360, 1]))
    
    # convert `radius` to iterable
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    
    if dxdy is None:
        dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
        dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
    else:
        dx, dy = dxdy
        
    nums_t = theta.size
    nums_r = radius.size
    
    # the meridional/zonal distance between center and every grid points
    dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    dist_lat = latlon2distance(lon, clat, lon, lat)
    dist_full = np.sqrt(dist_lon**2 + dist_lat**2)
    
    # give signs to dist_lon/lat (west of `clon` or south of `clat` would be negative)
    dist_lon[lon < clon] *= -1
    dist_lat[lat < clat] *= -1
    
    # set a box area to reduce the amount of computation
    max_r = radius.max()
    L_lon = int(max_r // dx) + 6
    L_lat = int(max_r // dy) + 6
    cix = np.unravel_index(dist_full.argmin(), dist_full.shape)  # (nearly) center index
    dist_lon_b = dist_lon[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
    dist_lat_b = dist_lat[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
    
    # reshape into interpolation form
    dist_lonlat_b = np.vstack((dist_lon_b.ravel(), dist_lat_b.ravel())).T
    
    # construct circular samples points and interpolate
    circle_pts = np.zeros((nums_t*nums_r, 2))
    circle_pts[:,0] = radius.repeat(nums_t) * np.cos(np.tile(theta, nums_r))
    circle_pts[:,1] = radius.repeat(nums_t) * np.sin(np.tile(theta, nums_r))
    interp_obj = FastGriddata(dist_lonlat_b, circle_pts)
    
    def inner(values):
        """
        Parameters:
        ----------
        values: 1-d array. see `circular_avg`
        
        Returns:
        -------
        circular mean
        """
        values_b = values[cix[0]-L_lat:cix[0]+L_lat, cix[1]-L_lon:cix[1]+L_lon]
        val_b = values_b.ravel()
        interp = interp_obj.interpolate(val_b)
        return interp.reshape(nums_r, nums_t).mean(axis=1)
    
    return inner


def rmw(lon, lat, ws, clon, clat, maxdist=550, dr=1, box=True, **kwargs):
    """
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
        The maximum search distance (km). Default is 550
    dr : scalar
        The radius (km) interval. Default is 1
    box : bool
        Create a box area that only calculate interpolation in this
        box. Default is True.
    **kwargs : 
        Keyword arguments for `circular_avg`
        
    Return
    ------
    radius : 1d array, shape = (n,)
        Radius, the coordinate of `axissym_ws`.
        The distance between i'th `axissym_ws` profile location
        and TC center is radius[i].
    axissym_ws : 1d array, shape = (n,)
        Axis-symmetric wind speed profile
    rmw : scalar
        Radius of maximum wind speed
    """
    if box:
        # find the index where lon[idx1,idx2] / lat[idx1,idx2] are the nearest to clon/clat
        _, idx2 = np.unravel_index(np.argmin(np.abs((lon - clon))), lon.shape)
        idx1, _ = np.unravel_index(np.argmin(np.abs((lat - clat))), lat.shape)

        # set a box area, and only calculate curcular average in this box area
        dlon = lon[0,1] - lon[0,0]
        dlat = lat[1,0] - lat[0,0]
        dy = dlat * 110.567
        L = (maxdist+10) // dy   # the half length of box edge

        bottom = (idx1 - L) if (idx1 - L >= 0) else 0
        up = (idx1 + L) if (idx1 + L <= lon.shape[0]) else -1
        left = (idx2 - L) if (idx2 - L >= 0) else 0
        right = (idx2 + L) if (idx2 + L <= lon.shape[1]) else -1
        bottom, up, left, right = int(bottom), int(up), int(left), int(right)

        lon_box = lon[bottom:up, left:right]
        lat_box = lat[bottom:up, left:right]
        ws_box = ws[bottom:up, left:right]
    else:
        lon_box = lon
        lat_box = lat
        ws_box = ws

    #radius = np.linspace(0, maxdist, maxdist)
    radius = np.arange(0, maxdist+dr, dr)
    axissym_ws = circular_avg(lon_box, lat_box, ws_box, clon, clat, radius)
    rmw = radius[np.nanargmax(axissym_ws)]
    return radius, axissym_ws, rmw


def axisymmetricity(lon, lat, var, radius, clon, clat, dxdy=None, integ='trapz', **kwargs):
    """
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
    **kwargs : 
        Keyword arguments for `circular_avg`
        
    Return:
    ------
    1d array with shape = (nr,). The axisymmetricity at the given radius.
    
    Reference
    ---------
    [1] Yoshiaki Miyamoto and Tetsuya Takemi: "A Transition Mechanism for the Spontaneous 
        Axisymmetric Intensification of Tropical Cyclones"
        J. Atmos. Sci, 70, 112-129
        https://doi.org/10.1175/JAS-D-11-0285.1
    """
    if dxdy is None:
        dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
        dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
        dxdy = (dx, dy)
        
    #if kwargs.get('theta'):
    #    _, _, dtheta = kwargs.get('theta')
    #else:
    #    dtheta = 2 * np.pi / 360   # the default `dtheta` in `circular_avg`
    theta = kwargs.pop('theta', np.arange(*np.deg2rad([0, 360, 1])))
    dtheta = theta[1] - theta[0]
        
    if integ == 'trapz':
        int_func = trapz
    elif integ == 'simps':
        int_func = simps
    else:
        raise ValueError(f'Unavailable `integ`: {integ}. It should be "trapz" or "simps".')
        
    kwargs.pop('return_interp', None)
    #interp = circular_avg(lon, lat, var, clon, clat, radius, dxdy, return_interp=True, **kwargs)
    interp = interp_circle(lon, lat, var, clon, clat, radius, theta, dxdy, coord='lonlat', **kwargs)
    ciravg = interp.mean(axis=1)    # (n_radius,)
    cirdev = interp - ciravg[:,np.newaxis]    # (n_radius, n_theta)
        
    res = ciravg**2 / (ciravg**2 + int_func(cirdev**2, dx=dtheta, axis=1) / (2*np.pi))   # (n_radius,)
    return res