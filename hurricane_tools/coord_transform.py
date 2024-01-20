import numpy as np

from . import config
from .distance import latlon2distance
from .pseudo_coord import lonlat2xy
from .interpolate import FastGriddata

try:
    from .fortran.f90xy2rt import xy2rt as f90xy2rt
    _HAS_IMPORT_F90XYRT = True
except ModuleNotFoundError:
    _HAS_IMPORT_F90XYRT = False
    
try:
    from .fortran.f90rt2xy import rt2xy as f90rt2xy
    _HAS_IMPORT_F90RT2XY = True
except ModuleNotFoundError:
    _HAS_IMPORT_F90RT2XY = False

try:
    from .fortran.f90interpz import mod_interpz as f90intp
    #from .fortran.f90interpz import find_level_1, find_level_n, calc_weights_1, calc_weights_n, interpz3d_1, interpz3d_n
    _HAS_IMPORT_F90INTERPZ = True
except ModuleNotFoundError:
    _HAS_IMPORT_F90INTERPZ = False


__all__ = [
    'XY2RT',
    'RT2XY',
    'Interpz3d'
]


def _convert_dtype(arr, dtype):
    if arr.dtype != np.dtype(dtype):
        arr = arr.astype(dtype, copy=False)
    return arr


class XY2RT:
    """
    Coordinate transformation from cartesian (x-y) to polar (radius-theta) coordinate.
    
    Note
    ----
    When the fortran routine is used, the performance will be further better if the 
    data types of all numpy array are `np.float32`.
    
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
        Comparison between intp='griddata' and 'fortran'
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
        lon = _convert_dtype(lon, 'float32')
        lat = _convert_dtype(lat, 'float32')
        
        if theta is None:
            theta = np.deg2rad(np.arange(360), dtype=np.float32)
            
        # convert `radius` to be iterable
        if isinstance(radius, (int, float, np.integer, np.floating)):
            radius = np.array([radius], dtype=np.float32)
        elif isinstance(radius, (list, tuple)):
            radius = np.array(radius, dtype=np.float32)
        elif isinstance(radius, np.ndarray):
            radius = _convert_dtype(radius, 'float32')
            
        if dxdy is None:
            dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
            dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
            dxdy = (dx, dy)
               
        # check interpolate option
        if intp == 'griddata':
            self._intp = 'griddata'
        elif intp == 'fortran':
            self._intp = 'fortran'
        elif intp is None:
            if config.INTERP_OPTION == 'griddata':
                self._intp = 'griddata'
            elif config.INTERP_OPTION == 'fortran':
                self._intp = 'fortran'
            else:
                raise ValueError("`config.INTERP_OPTION` must be 'griddata' or 'fortran'.")
        else:
            raise ValueError("`intp` must be None, 'griddata' or 'fortran'.")
            
        if (self._intp == 'fortran') and (not _HAS_IMPORT_F90XYRT):
            raise ModuleNotFoundError("The Fortran modules haven't been compiled.")
            
        if self._intp == 'griddata':
            self._intp_func = self._prepare_xy2rt_griddata(lon, lat, clon, clat, radius, theta, dxdy)
        elif self._intp == 'fortran':
            self._intp_func = self._prepare_xy2rt_f90(lon, lat, clon, clat, radius, theta, dxdy)
            
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
        if self._intp == 'fortran':
            var = _convert_dtype(var, 'float32')
        return self._intp_func(var)
        
    def _prepare_xy2rt_griddata(self, lon, lat, clon, clat, radius, theta, dxdy):
        """return a function to perform interpolation by griddata method"""
        dx, dy = dxdy
        ny, nx = lon.shape
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
        lowery = cix[0] - L_lat if (cix[0] - L_lat > 0) else 0
        uppery = cix[0] + L_lat if (cix[0] + L_lat < ny) else ny
        lowerx = cix[1] - L_lon if (cix[1] - L_lon > 0) else 0
        upperx = cix[1] + L_lon if (cix[1] + L_lon > nx) else nx
        
        box_slice = (slice(lowery, uppery), slice(lowerx, upperx))
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
    
    Note
    ----
    When the fortran routine is used, the performance will be further better if the 
    data types of all numpy array are `np.float32`.
    
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
        lon = _convert_dtype(lon, 'float32')
        lat = _convert_dtype(lat, 'float32')
        radius = _convert_dtype(radius, 'float32')
        theta = _convert_dtype(theta, 'float32')
        
        theta = np.mod(theta, 2*np.pi)
        
        # check interpolate option
        if intp == 'griddata':
            self._intp = 'griddata'
        elif intp == 'fortran':
            self._intp = 'fortran'
        elif intp is None:
            if config.INTERP_OPTION == 'griddata':
                self._intp = 'griddata'
            elif config.INTERP_OPTION == 'fortran':
                self._intp = 'fortran'
            else:
                raise ValueError("`config.INTERP_OPTION` must be 'griddata' or 'fortran'.")
        else:
            raise ValueError("`intp` must be None, 'griddata' or 'fortran'.")
            
        if (self._intp == 'fortran') and (not _HAS_IMPORT_F90XYRT):
            raise ModuleNotFoundError("The Fortran modules haven't been compiled.")
            
        if self._intp == 'griddata':
            self._intp_func = self._prepare_rt2xy_griddata(lon, lat, clon, clat, radius, theta)
        elif self._intp == 'fortran':
            self._intp_func = self._prepare_rt2xy_f90(lon, lat, clon, clat, radius, theta)
    
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
        if self._intp == 'fortran':
            var = _convert_dtype(var, 'float32')
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
    
    Note
    ----
    The performance will be better if the data types of numpy arrays are all `np.float32`
    or all `np.float64`.
    If the types of arrays are not all identical, or the types are not `np.float32` or
    `np.float64`, implicitly type conversion is needed and the performance will be affected.
    
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
        if not _HAS_IMPORT_F90INTERPZ:
            raise ModuleNotFoundError("The fortran modules haven't been compiled.")
        
        self._dtype = zdata.dtype
        self.zdata = zdata
        self.level = np.array(level, dtype=self._dtype)
        self.missing_value = missing_value
        self._nz = zdata.shape[0]
        
        # the used functions depend on `levels` (single or multiple) and float precision (single or double)
        if isinstance(level, (int, float, np.integer, np.floating)):
            postfix1 = '1'
        elif isinstance(level, (tuple, list, np.ndarray)):
            postfix1 = 'n'
        else:
            raise ValueError(f"Unavailable `level` type : {type(level)}")
        
        if self._dtype == 'float32':
            self._is_convert_dtype = False
            postfix2 = 'r32'
        elif self._dtype == 'float64':
            self._is_convert_dtype = False
            postfix2 = 'r64'
        else:
            self._is_convert_dtype = True
            self.zdata = self.zdata.astype(np.float32)
            postdix2 = 'r32'
            
        postfix = f'{postfix1}_{postfix2}'
        find_level_func = getattr(f90intp, f'find_level_{postfix}')
        find_weight_func = getattr(f90intp, f'calc_weights_{postfix}')
        self._interpz3d_func = getattr(f90intp, f'interpz3d_{postfix}')
        
        # find level index and weights
        zdata_f = np.asfortranarray(zdata.T)
        lev_idx = find_level_func(zdata_f, level)    # (nx, ny, nlev)
        weight = find_weight_func(zdata_f, level, lev_idx)   # (nx, ny, nlev)
        
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
        dtype = self._dtype
        nz = self._nz
        lev_idx = self._lev_idx
        weight = self._weight
        
        interpz3d_func = self._interpz3d_func
        missing_value = self.missing_value
        
        if len(var) == 1:
            if self._is_convert_dtype:
                v = _convert_dtype(var[0], dtype)
            else:
                v = var[0]
            
            if v.shape != self.zdata.shape:
                raise ValueError("The shape of input array should be equal to `zdata`.")
            
            v_f = np.asfortranarray(v.T)
            v_interp = interpz3d_func(v_f, lev_idx, weight)   # (nx, ny)
            v_interp = np.where((lev_idx == 0) | (lev_idx == nz), missing_value, v_interp)
            return v_interp.T
        
        else:
            var_interp = []
        
            for v in var:
                if self._is_convert_dtype:
                    v = _convert_dtype(v, dtype)
                
                if v.shape != self.zdata.shape:
                    raise ValueError("The shape of input array should be equal to `zdata`.")
                    
                v_f = np.asfortranarray(v.T)
                v_interp = interpz3d_func(v_f, lev_idx, weight)   # (nx, ny, nlev)
                v_interp = np.where((lev_idx == 0) | (lev_idx == nz), missing_value, v_interp)
                var_interp.append(v_interp.T)
            
            return var_interp