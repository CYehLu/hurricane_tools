import numpy as np
from scipy.interpolate import griddata
from scipy.integrate import trapz, simps

from .distance import latlon2distance


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
    If `return_interp` is False (default option):
        Circular average result, 1-d array with size=len(radius).
    If `return_interp` is True:
        See `return_interp` argument interpretation above.
        
    NOTE:
    ----
    If this function is used repeatly, and all parameters remain the same except `values`, it is 
    worth to use `circular_avg_closure` instead. `circular_avg_closure` returns a closure function
    which only use `values` as its parameter.
    """
    # check kwargs, the remainning part of kwargs would only contain the arguments of `griddata`
    if kwargs.get('theta'):
        theta_start, theta_end, dtheta = kwargs.get('theta')
        thetas = np.arange(theta_start, theta_end, dtheta)
        kwargs.pop('theta', None)
    else:
        theta_start, theta_end, dtheta = np.deg2rad([0, 360, 1])
        thetas = np.arange(theta_start, theta_end, dtheta)
    nums_t = thetas.size
    
    if kwargs.get('return_interp'):
        return_flag = True
    else:
        return_flag = False
    kwargs.pop('return_interp', None)
        
    # check `dxdy`, if None then determine `dx` and `dy` based on `lon` and `lat`
    if dxdy is None:
        dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
        dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
    else:
        dx, dy = dxdy
    
    # convert `radius` to iterable
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    nums_r = radius.size
    
    # the meridional/zonal distance between center and every grid points
    dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    dist_lat = latlon2distance(lon, clat, lon, lat)    # (ny, nx)
    dist_full = np.sqrt(dist_lon**2 + dist_lat**2)    # (ny, nx)
    
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
    circle_pts[:,0] = radius.repeat(nums_t) * np.cos(np.tile(thetas, nums_r))
    circle_pts[:,1] = radius.repeat(nums_t) * np.sin(np.tile(thetas, nums_r))
    interp = griddata(dist_lonlat_b, val_b, circle_pts, **kwargs)
    
    if return_flag:
        return interp.reshape(nums_r, nums_t)
    else:
        res = interp.reshape(nums_r, nums_t).mean(axis=1)
        return res


def circular_avg_closure(lon, lat, clon, clat, radius, dxdy=None, **kwargs):
    """
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
    """
    # check kwargs, the remainning part of kwargs would only contain the arguments of `griddata`
    if kwargs.get('theta'):
        theta_start, theta_end, dtheta = kwargs.get('theta')
        thetas = np.arange(theta_start, theta_end, dtheta)
        kwargs.pop('theta', None)
    else:
        theta_start, theta_end, dtheta = np.deg2rad([0, 360, 1])
        thetas = np.arange(theta_start, theta_end, dtheta)
    nums_t = thetas.size
    
    if kwargs.get('return_interp'):
        return_flag = True
    else:
        return_flag = False
    kwargs.pop('return_interp', None)
        
    # check `dxdy`, if None then determine `dx` and `dy` based on `lon` and `lat`
    if dxdy is None:
        dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
        dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
    else:
        dx, dy = dxdy
    
    # convert `radius` to iterable
    if isinstance(radius, (int, float)):
        radius = np.array([radius])
    elif isinstance(radius, list):
        radius = np.array(radius)
    nums_r = radius.size
    
    # the meridional/zonal distance between center and every grid points
    dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
    dist_lat = latlon2distance(lon, clat, lon, lat)    # (ny, nx)
    dist_full = np.sqrt(dist_lon**2 + dist_lat**2)    # (ny, nx)
    
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
    
    # construct circular sample points
    circle_pts = np.zeros((nums_t*nums_r, 2))
    circle_pts[:,0] = radius.repeat(nums_t) * np.cos(np.tile(thetas, nums_r))
    circle_pts[:,1] = radius.repeat(nums_t) * np.sin(np.tile(thetas, nums_r))
    
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
        interp = griddata(dist_lonlat_b, val_b, circle_pts, **kwargs)
        
        if return_flag:
            return interp.reshape(nums_r, nums_t)
        else:
            res = interp.reshape(nums_r, nums_t).mean(axis=1)
            return res
    
    return inner


def rmw(lon, lat, ws, clon, clat, maxdist=550, box=True, **kwargs):
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
        
    # 2020/07/01: update this part because of the update of `circular_avg`
    """
    # prepare to circular average
    maxlat = maxdist / 110.567   # 1 lat degree = 110.567 km
    lonlat = np.hstack([lon_box.reshape(-1, 1), lat_box.reshape(-1, 1)])
    ws_1d = ws_box.ravel()
    r = np.linspace(0, maxlat, maxdist) 
    
    # find the axissymmetric wind profile and RMW
    axissym_ws = circular_avg(lonlat, ws_1d, clon, clat, r, **kwargs)
    rmw = r[np.nanargmax(axissym_ws)] * 110.567
    dist_coord = r * 110.567
    
    return dist_coord, axissym_ws, rmw
    """
    
    radius = np.linspace(0, maxdist, maxdist)
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
    """
    if dxdy is None:
        dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
        dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
        dxdy = (dx, dy)
        
    if kwargs.get('theta'):
        _, _, dtheta = kwargs.get('theta')
    else:
        dtheta = 2 * np.pi / 360   # the default `dtheta` in `circular_avg`
        
    if integ == 'trapz':
        int_func = trapz
    elif integ == 'simps':
        int_func = simps
    else:
        raise ValueError(f'Unavailable `integ`: {integ}. It should be "trapz" or "simps".')
        
    kwargs.pop('return_interp', None)
    interp = circular_avg(lon, lat, var, clon, clat, radius, dxdy, return_interp=True, **kwargs)
    ciravg = interp.mean(axis=1)    # (n_radius,)
    cirdev = interp - ciravg[:,np.newaxis]    # (n_radius, n_theta)
        
    res = ciravg**2 / (ciravg**2 + int_func(cirdev**2, dx=dtheta, axis=1) / (2*np.pi))   # (n_radius,)
    return res