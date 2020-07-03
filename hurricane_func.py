import pandas as pd
import numpy as np
from scipy.optimize import root
from scipy.interpolate import griddata
from scipy.integrate import trapz, simps


def find_centers_nearest_time(filename, target_dates, formats=None):
    """
    Find the hurricane center data at specific times from the file
    which contains center information based on WRF rsl.error.0000.
    
    Parameters:
    ----------
    filename: 
        str, the file name of WRF hurricane centers file.
        The file can be obtained by exexute 'get_wrf_ty_centers'.
    target_dates: 
        list of pandas datetime objects or list of str.
        The nearest datetime and its data of the center data in the 'filename'
        file would be find out.
    formats:
        str. The string format of target_dates if they are str. Option.
        
    Return:
    ------
    DataFrame
    """
    # read center file
    with open(filename) as cf:
        contents = cf.readlines()
        
    # only keep the time, lat and lon information
    contents = [line.split()[1:4] for line in contents]
    
    # parse contents
    df = pd.DataFrame(contents)
    df = df.set_index(0)
    df.index.name = None
    df.columns = ['lat', 'lon']
    df.index = pd.to_datetime(df.index, format='%Y-%m-%d_%H:%M:%S')
    df = df.applymap(lambda s: float(s))
    
    # convert to list of pd.datetime obj if "target_dates" are strings
    if isinstance(target_dates[0], str):
        target_dates = [pd.to_datetime(dt) for dt in target_dates]
        
    # find nearest dates
    idx = []
    for dt in target_dates:
        idx.append(df.index.get_loc(dt, method='nearest'))
        
    df2 = pd.DataFrame(
        data=df.iloc[idx,:].values,
        index=target_dates,
        columns=['lat', 'lon']
    )

    return df2


def latlon2distance(lon1, lat1, lon2, lat2):
    """calculate the distance (km) of two positions"""
    # approximate radius of earth in km
    R = 6373.0
    
    # convert to radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    return R * c


def find_lonlat_with_distance(clon, clat, distance, xy):
    """
    Find the longtitude/latitude which the horizontal or vertical distance
    to the (clon, clat) is equal to 'distance'.
    
    Parameters:
    ----------
    clon, clat: scaler. The reference coordinate.
    distance: scaler. The distance (km) between target coordinate and (clon, clat).
    xy: str, 'x' or 'y'. 
        If x, it will find the coordinate which distance along latitude line is
        equal to 'distance' (The horizontal distance).
        If y, it will find the coordinate which distance along longtitude line is
        equal to 'distance' (The vertical distance).
        
    Return:
    ------
    scaler, the difference of lon/lat degree
    """
    if xy == 'x':
        f = lambda ll: latlon2distance(clon, clat, ll, clat) - distance
        return np.abs(root(f, clon).x[0] - clon)
    elif xy == 'y':
        f = lambda ll: latlon2distance(clon, clat, clon, ll) - distance
        return np.abs(root(f, clat).x[0] - clat)
    
    
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
        The radius of circle. Unit is km.
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
    circular average result, 1-d array with size=len(radius).
        
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
    
    
def ty_center_mslp(lon, lat, slp):
    """
    Find typhoon center by minimum sea level pressure grid (no weighted).
    
    Parameter
    ---------
    lon, lat : 2d numpy array, shape = (ny, nx)
        longtitude / latitude
    slp : 2d numpy array, shape = (ny, nx)
        sea level pressure
        
    Return
    ------
    tuple, (center_lon, center_lat)
    """
    i, j = np.unravel_index(np.argmin(slp), slp.shape)
    clon, clat = lon[i,j], lat[i,j]
    return clon, clat
    
    
def weighted_ty_center(lon, lat, center_lon, center_lat, var, latlon_range=1):
    """
    Giving first guess typhoon center longtitude and latitude, calculate the 
    new center by weighted mothod.
    
    Parameter:
    ---------
    lon, lat: 
        2-d numpy ndarray. Their shape should be equal.
    center_lon, center_lat: 
        float, the first guess typhoon center coordinate.
    var:
        2-d numpy ndarray. Used for determination the center of typhoon, usually
        is pressure.
        Its shape should equal to lon and lat.
    latlon_range: 
        scaler, the degree of lat/lon which will generate a box to calculate the
        weights.
        
    Return:
    ------
    tuple, (new_lon, new_lat). 
    """
    # find the box
    bool_lat = (lat >= center_lat - latlon_range) & (lat <= center_lat + latlon_range)
    bool_lon = (lon >= center_lon - latlon_range) & (lon <= center_lon + latlon_range)
    bool_lonlat = bool_lat & bool_lon
    where = np.where(bool_lonlat)
    left_grid = np.min(where[0])
    right_grid = np.max(where[0])
    bottom_grid = np.min(where[1])
    top_grid = np.max(where[1])
    
    # restrict the lat, lon and var in the box range
    box_lat = lat[left_grid:right_grid,bottom_grid:top_grid]
    box_lon = lon[left_grid:right_grid,bottom_grid:top_grid]
    box_var = var[left_grid:right_grid,bottom_grid:top_grid]
    
    mean_box_var = np.mean(box_var)
    
    # only the positive deviation values can be kept
    weights = np.abs(mean_box_var - box_var) * np.heaviside(mean_box_var - box_var, 0)
    new_center_lon = np.sum(weights * box_lon) / np.sum(weights)
    new_center_lat = np.sum(weights * box_lat) / np.sum(weights)
    
    return (new_center_lon, new_center_lat)


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


def axissymmetricity(lon, lat, var, radius, clon, clat, dxdy=None, integ='trapz', **kwargs):
    """
    var: shape=(ny, nx)
    radius: shape=(nr,)
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
