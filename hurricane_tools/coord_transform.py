import numpy as np

from . import config
from .distance import latlon2distance
from .pseudo_coord import lonlat2xy
from .interpolate import FastGriddata

try:
    from .fortran.f90xy2rt import xy2rt as f90xy2rt
    _has_import_f90xy2rt = True
except ModuleNotFoundError:
    _has_import_f90xy2rt = False
    
try:
    from .fortran.f90rt2xy import rt2xy as f90rt2xy
    _has_import_f90rt2xy = True
except ModuleNotFoundError:
    _has_import_f90rt2xy = False

try:
    from .fortran.f90interpz import find_level_1, find_level_n, calc_weights_1, calc_weights_n, interpz3d_1, interpz3d_n
    _has_import_f90interpz = True
except ModuleNotFoundError:
    _has_import_f90interpz = False


__all__ = [
    'XY2RT',
    'RT2XY',
    'Interpz3D'
]


class XY2RT:
    """
    Coordinate transformation from cartesian (x-y) to polar (radius-theta) coordinate.
    
    Example
    -------
    >>> from hurricane_tools import coord_transform
    >>> lon, lat, clon, clat, var = ...   # get data
    >>> var.shape       # (ny, nx)
    (200, 200)
    >>> radius = np.arange(0, 100, 2)
    >>> xy2rt = coord_transform.XY2RT(lon, lat, clon, clat, radius)
    >>> var_rt = xy2rt(var)
    >>> var_rt.shape    # (nradius, ntheta)
    (50, 360)
    """
    
    def __init__(self, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None):
        """
        Initialization.
        
        Parameters
        ----------
        lon, lat: 2-d array, shape = (ny, nx)
            The longitude/latitude coordinates
        clon, clat: scalar
            Center coordinate
        radius: scalar or 1-d array-like, shape = (nradius,)
            The radius to be interpolated. Unit is km.
        theta: 1-d array, shape = (ntheta,). Optional
            The azimuth to be interpolated. Unit is radians.
            Default is np.deg2rad(np.arange(360))
            It should be "ascent" order (clockwise), and the values should be greater than 0.
        dxdy: Tuple(scalar, scalar). Optional
            Spatial resolution, (dx, dy). Unit is km.
            Default is None, and it would be automatically derived based on `lon` and `lat`.
        intp: str, 'griddata' or 'fortran'. Optional
            Choose the interpolation method, the first one is based on `hurricane_tools.interpolate.FastGriddata`
            and the second one is based on the fortran subroutine.
            The difference between these two methods can be seen in `Note` below.
            If None (default), it would use the setting in `hurricane_tools.config.INTERP_OPTION`.
            
        Note
        ----
        griddata: 
            Slower, need more data points around the interpolation point.
            Slightly more accurate.
            Only depends on python/numpy/scipy.
        fortran:
            Faster, only need 4 data points around the interpolation point.
            Slightly less accurate.
            Need to compile fortran code before using it (I only test it on linux system with
            gfortran compiler).
        """
        
        if theta is None:
            theta = np.deg2rad(np.arange(360))
            
        # convert `radius` to be iterable
        if isinstance(radius, (int, float)):
            radius = np.array([radius])
        elif isinstance(radius, (list, tuple)):
            radius = np.array(radius)
            
        if dxdy is None:
            dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
            dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
            dxdy = (dx, dy)
        
        # check interpolate option
        if intp is None:
            if config.INTERP_OPTION == 'griddata':
                self._intp_func = self._prepare_xy2rt_griddata(lon, lat, clon, clat, radius, theta, dxdy)
                
            elif config.INTERP_OPTION == 'fortran':
                if _has_import_f90xy2rt:
                    self._intp_func = self._prepare_xy2rt_f90(lon, lat, clon, clat, radius, theta, dxdy)
                else:
                    raise ModuleNotFoundError("The fortran modules haven't been compiled.")
                    
            else:
                raise ValueError("`config.INTERP_OPTION` must be 'griddata' or 'fortran'.")
            
        elif intp == 'griddata':
            self._intp_func = self._prepare_xy2rt_griddata(lon, lat, clon, clat, radius, theta, dxdy)
        
        elif intp == 'fortran':
            if _has_import_f90xy2rt:
                self._intp_func = self._prepare_xy2rt_f90(lon, lat, clon, clat, radius, theta, dxdy)
            else:
                raise ModuleNotFoundError("The fortran modules haven't been compiled.")
                
        else:
            raise ValueError("`intp` must be None, 'griddata' or 'fortran'.")
    
    def __call__(self, var):
        """
        Interpolate data from cartesian (x-y) to polar (radius-theta) coordinate.
        
        Parameter
        ---------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Interpolating result, shape = (nradius, ntheta)
        """
        return self._intp_func(var)
        
    def _prepare_xy2rt_griddata(self, lon, lat, clon, clat, radius, theta, dxdy):
        """return a function to perform interpolation by griddata method"""
        dx, dy = dxdy
        nrad = radius.size
        nthe = theta.size

        # the meridional/zonal distance between center and grid points
        dist_lon = latlon2distance(clon, lat, lon, lat)    # (ny, nx)
        dist_lat = latlon2distance(lon, clat, lon, lat)
        dist_full = np.sqrt(dist_lon**2 + dist_lat**2)

        # signed distance
        dist_lon[lon < clon] *= -1
        dist_lat[lat < clat] *= -1

        # set a box area to reduce the amount of computation
        max_rad = radius.max()
        L_lon = int(max_rad // dx) + 6
        L_lat = int(max_rad // dy) + 6
        cix = np.unravel_index(dist_full.argmin(), dist_full.shape)  # (nearly) center index
        box_slice = (slice(cix[0]-L_lat, cix[0]+L_lat), slice(cix[1]-L_lon, cix[1]+L_lon))
        dist_lon_b = dist_lon[box_slice]
        dist_lat_b = dist_lat[box_slice]

        # reshape to interpolation form
        dist_lonlat_b = np.vstack((dist_lon_b.ravel(), dist_lat_b.ravel())).T

        # construct circular samples points and interpolate
        circle_pts = np.zeros((nthe*nrad, 2))
        circle_pts[:,0] = radius.repeat(nthe) * np.cos(np.tile(theta, nrad))
        circle_pts[:,1] = radius.repeat(nthe) * np.sin(np.tile(theta, nrad))

        # FastGriddata instance
        fgd_obj = FastGriddata(dist_lonlat_b, circle_pts)

        def _inner(values):
            values_b = values[box_slice]
            val_b = values_b.ravel()
            interp = fgd_obj.interpolate(val_b)
            return interp.reshape(nrad, nthe)

        return _inner
    
    def _prepare_xy2rt_f90(self, lon, lat, clon, clat, radius, theta, dxdy):
        """return a function to perform interpolation by fortran method"""
        dx, dy = dxdy
        f90xy2rt.undef = np.nan
        
        grididx, weights = f90xy2rt.get_grididx_and_weights(
            np.asfortranarray(lon), 
            np.asfortranarray(lat), 
            clon, clat, 
            radius, theta,
            dx, dy
        )
        
        def _inner(values):
            return f90xy2rt.interp_xy2rt(grididx, weights, np.asfortranarray(values))
        
        return _inner
    
    
class RT2XY:
    """
    Coordinate transformation from polar (radius-theta) to cartesian (x-y) coordinate.
    
    Example
    -------
    >>> from hurricane_tools import coord_transform
    >>> lon, lat, clon, clat, var = ...   # get data
    >>> var.shape       # (ny, nx)
    (200, 200)
    >>> radius = np.arange(0, 100, 2)
    >>> theta = np.deg2rad(np.arange(360))
    >>> xy2rt = coord_transform.XY2RT(lon, lat, clon, clat, radius, theta)
    >>> var_rt = xy2rt(var)   # shape = (nradius, ntheta)
    >>> rt2xy = coord_transform.RT2XY(lon, lat, clon, clat, radius, theta)
    >>> var_rt2xy = rt2xy(var_rt)
    >>> var_rt2xy.shape    # shape = (ny, nx)
    (200, 200)
    """
    
    def __init__(self, lon, lat, clon, clat, radius, theta, intp=None):
        """
        Initialization.
        
        Parameters
        ----------
        lon, lat: 2-d array, shape = (ny, nx)
            The longitude/latitude coordinates to be interpolated.
        clon, clat: scalar
            Center coordinate
        radius: scalar or 1-d array-like, shape = (nradius,)
            The radius coordinate of data. Unit is km.
        theta: 1-d array, shape = (ntheta,)
            The azimuth coordinate of data. Unit is radians.
            It should be "ascent" order (clockwise), and the values should be greater than 0.
        intp: str, 'griddata' or 'fortran'. Optional
            Choose the interpolation method, the first one is based on `interpolate.FastGriddata`
            and the second one is based on the fortran subroutine.
            The difference between these two methods can be seen in `Note` in `XY2RT`.
            If None (default), it would use the setting in `config.INTERP_OPTION`.
        """
        theta = np.mod(theta, 2*np.pi)
        
        if intp is None:
            if config.INTERP_OPTION == 'griddata':
                self._intp_func = self._prepare_rt2xy_griddata(lon, lat, clon, clat, radius, theta)
                
            elif config.INTERP_OPTION == 'fortran':
                if _has_import_f90rt2xy:
                    self._intp_func = self._prepare_rt2xy_f90(lon, lat, clon, clat, radius, theta)
                else:
                    raise ModuleNotFoundError("The fortran modules haven't been compiled.")
                    
            else:
                raise ValueError("`config.INTERP_OPTION` must be 'griddata' or 'fortran'.")
            
        elif intp == 'griddata':
            self._intp_func = self._prepare_rt2xy_griddata(lon, lat, clon, clat, radius, theta)
        
        elif intp == 'fortran':
            if _has_import_f90rt2xy:
                self._intp_func = self._prepare_rt2xy_f90(lon, lat, clon, clat, radius, theta)
            else:
                raise ModuleNotFoundError("The fortran modules haven't been compiled.")
                
        else:
            raise ValueError("`intp` must be None, 'griddata' or 'fortran'.")
    
    def __call__(self, var):
        """
        Interpolate data from polar (radius-theta) to cartesian (x-y) coordinate.
        
        Parameter
        ---------
        var: 2-d array, shape = (nradius, ntheta)
        
        Return
        ------
        Interpolating result, shape = (ny, nx)
        """
        return self._intp_func(var)
    
    def _prepare_rt2xy_griddata(self, lon, lat, clon, clat, radius, theta):
        """return a function to perform interpolation by griddata method"""
        X, Y = lonlat2xy(lon, lat, clon, clat)
        
        t_mesh, r_mesh = np.meshgrid(theta, radius)    # shape = (nrad, nthe)
        polar_coord = np.vstack((r_mesh.ravel(), t_mesh.ravel())).T
        
        interp_coord_r = np.sqrt(X**2 + Y**2)
        interp_coord_t = np.arctan2(Y, X)
        interp_coord_t[interp_coord_t < 0] += 2*np.pi
        interp_coord = np.vstack((interp_coord_r.ravel(), interp_coord_t.ravel())).T
        
        fgd_obj = FastGriddata(polar_coord, interp_coord)
        
        def _inner(values):
            return fgd_obj.interpolate(values.ravel()).reshape(X.shape)
        
        return _inner
    
    def _prepare_rt2xy_f90(self, lon, lat, clon, clat, radius, theta):
        """return a function to perform interpolation by fortran method"""
        f90rt2xy.undef = np.nan
        
        radthe_idx, radthe = f90rt2xy.get_radtheidx_radthe(
            np.asfortranarray(lon), 
            np.asfortranarray(lat),
            clon, clat,
            radius, theta
        )
        
        def _inner(values): 
            return f90rt2xy.interp_rt2xy(radius, theta, radthe_idx, radthe, np.asfortranarray(values))
        
        return _inner


class Interpz3d:
    """
    Interpolating variables on specified vertical coordinate.
    
    Example
    -------
    >>> u, v, phi, temp, pres = get_data()          # a fake function to get data
    >>> interp_obj = Interpz3d(pres, [900, 850])    # interpolate variables on 900 and 850 hPa
    >>> u850, v850, phi850, temp850 = interp_obj.interp(u, v, phi, temp)   
    >>> u850.shape     # (2, ny, nx)
    """
    
    def __init__(self, zdata, level, missing_value=np.nan):
        """
        Initialize with pressure and levels.
        
        Parameter
        ---------
        zdata : 3-d array, shape = (nz, ny, nx)
            The variable of new vertical coordinate. e.g pressure or height
        level : scalar, or 1-d array-like with shape = (nlev,)
            The desired vertical levels. 
            For example, `Interpz3d(pressure, [900, 850])` indicates that the variables would be 
            interpolated on 900 and 850 hPa vertical levels.
            Note that extrapolation is not allowed. If some values in `level` are out of the range
            of `zdata`, the interpolated value will be assigned to `missing_value`.
        missing_value : scalar, optional
            Assign missing value. Default is `np.nan`
        """     
        if not _has_import_f90interpz:
            raise ModuleNotFoundError("The fortran modules haven't been compiled.")
        
        self.zdata = zdata
        self.level = np.array(level)
        self.missing_value = missing_value
        
        if isinstance(level, (int, float, np.integer, np.floating)):
            find_level_func = find_level_1
            find_weight_func = calc_weights_1
            self._interpz3d_func = interpz3d_1
        elif isinstance(level, (tuple, list, np.ndarray)):
            find_level_func = find_level_n
            find_weight_func = calc_weights_n
            self._interpz3d_func = interpz3d_n
        else:
            raise ValueError(f"Unavailable `level` type : {type(level)}")
            
        # find level index and weight
        zdata_f = np.asfortranarray(zdata.T)
        lev_idx = find_level_func(zdata_f, level)    # (nx, ny, nlev)
        weight = find_weight_func(zdata_f, level, lev_idx)    # (nx, ny, nlev)
        
        self._zdata_f = zdata_f
        self._lev_idx = lev_idx
        self._weight = weight
    
    def interp(self, *var):
        """
        Interpolate variable from the original vertical coordinate to `zdata` coordinate.
        
        Parameter
        ---------
        *var : 3-d array, shape = (nz, ny, nx)
            Any number of variables to be interpolated. Their shapes should be equal to `zdata`.
            
        Return
        ------
        If len(var) == 1:
            Interpolated variable. 3-d array, shape = (nlev, ny, nx)
        If len(var) > 1:
            List of interpolated variables. List[Array_3d], shape of each array = (nlev, ny, nx)
        """
        nz = self.zdata.shape
        lev_idx = self._lev_idx
        weight = self._weight
        
        interpz3d_func = self._interpz3d_func
        missing_value = self.missing_value
        
        if len(var) == 1:
            if var[0].shape != self.zdata.shape:
                raise ValueError("The shape of input array should be equal to `zdata`.")
            
            v_f = np.asfortranarray(var[0].T)
            v_interp = interpz3d_func(v_f, lev_idx, weight)   # (nx, ny)
            v_interp = np.where((lev_idx == 0) | (lev_idx == nz), missing_value, v_interp)
            return v_interp.T
        
        else:
            var_interp = []
        
            for v in var:
                if v.shape != self.zdata.shape:
                    raise ValueError("The shape of input array should be equal to `zdata`.")
                    
                v_f = np.asfortranarray(v.T)
                v_interp = interpz3d_func(v_f, lev_idx, weight)   # (nx, ny, nlev)
                v_interp = np.where((lev_idx == 0) | (lev_idx == nz), missing_value, v_interp)
                var_interp.append(v_interp.T)
            
            return var_interp