import warnings

import numpy as np
from scipy.integrate import trapz, simps

from .distance import latlon2distance
from .coord_transform import XY2RT, RT2XY


__all__ = [
    'CircularAvg',
    'CircularFourier',
    'rmw',
    'Axisymmetricity',
    'Rotate'
]


class CircularAvg(XY2RT):
    """
    Compute the azithmual mean field (axisymmetric field).
    """
    
    def __init__(self, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None):
        """
        Initialization.
        
        Parameter
        ---------
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
            See `hurricane_tools.coord_transform.XY2RT` for more details.
        """
        # add `self._intp_func()`
        super().__init__(lon, lat, clon, clat, radius, theta, dxdy, intp)
        
    def __call__(self, values):
        """
        Calculate the azimuthal mean.
        
        Parameter
        ---------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Azimuthal mean, shape = (nradius,)
        """
        return self._intp_func(values).mean(axis=1)
        
        
class CircularFourier:
    def __init__(self, lon, lat, clon, clat, radius, dxdy=None, intp=None):
        """
        Initialization.
        
        Parameter
        ---------
        lon, lat : 2-d array, shape = (ny, nx)
            Longtitude / Latitude
        clon, clat : scalar
            Central longtitude / latitude
        radius : 1-d array, shape = (nradius,)
            The radius of circles. Unit is km.
        dxdy : Tuple(scalar, scalar). Optional
            Spatial resolution, (dx, dy). Unit is km.
            Default is None, and it would be automatically derived based on `lon` and `lat`.
        intp: str, 'griddata' or 'fortran'. Optional
            Choose the interpolation method, the first one is based on `hurricane_tools.interpolate.FastGriddata`
            and the second one is based on the fortran subroutine.
            See `hurricane_tools.coord_transform.XY2RT` for more details.
        """
        self._dtheta = np.deg2rad(1)
        theta = np.deg2rad(np.arange(360))
        
        self._xy2rt_obj = XY2RT(lon, lat, clon, clat, radius, theta, dxdy, intp)
        self._rt2xy_obj = RT2XY(lon, lat, clon, clat, radius, theta, intp)
        
    def circular_fft(self, var):
        """
        Perform FFT on circles.
        
        Parameter
        ---------
        var : 2-d array, shape = (ny, nx)
            Variable to perform FFT
            
        Return
        ------
        k : 1-d array, shape = (m,)
            Wave numbers
        var_hat : 2-d array, shape = (nradius, m)
            FFT result
        """
        var_rt = self._xy2rt_obj(var)
        var_hat = np.fft.fft(var_rt, axis=1)    # (nradius, ntheta/m)
        k = 2 * np.pi * np.fft.fftfreq(var_rt.shape[1], d=self._dtheta)
        return k, var_hat
    
    def decomp(self, var, ns):
        """
        Decompute variable and extract the wave number `ns` components respectively.
        
        Parameter
        ---------
        var : 2-d array, shape = (ny, nx)
            Variable to decompute
        ns : List[int]
            The list of extracting wave numbers.
            
        Return
        ------
        result : List[ndarray], ndarray shape = (ny, nx)
            result[i] is the wave number `ns[i]` component of `var`.
        """           
        # perform fft on every radius
        k, var_hat = self.circular_fft(var)
        
        res = []
        for n in ns:
            # filter out all wave number except `n`
            var_hat_tmp = np.copy(var_hat)
            var_hat_tmp[:, ~(np.abs(k) == n)] = 0
            var_n = np.fft.ifft(var_hat_tmp, axis=1).real
            var_n_xy = self._rt2xy_obj(var_n)
            res.append(var_n_xy)
            
        return res


def rmw(lon, lat, ws, clon, clat, maxdist=None, mindist=None, dr=None, dxdy=None, intp=None):
    """
    Find the radius of maximum wind.

    Parameters
    ----------
    lon, lat: 2-d array, shape = (ny, nx)
        The longitude/latitude coordinates
    ws: 2-d array, shape = (ny, nx)
        Wind speed
    clon, clat: scalar
        Center coordinate
    maxdist: scalar. Optional
        The maximum search distance. Unit is km.
    mindist: scalar. Optional
        The minimum search distance. Unit is km.
    dr: scalar. Optional
        The interpolated radius interval.
        If None, it will be `max(dxdy)`.
    dxdy: Tuple(scalar, scalar). Optional
        Spatial resolution, (dx, dy). 
        Default is None, and it would be automatically derived based on `lon` and `lat`.
    intp: str, 'griddata' or 'fortran'. Optional
        Choose the interpolation method, the first one is based on `hurricane_tools.interpolate.FastGriddata`
        and the second one is based on the fortran subroutine.
        See `hurricane_tools.coord_transform.XY2RT` for more details.
    """
    if dxdy is None:
        dx = latlon2distance(lon[:,1:], lat[:,1:], lon[:,:-1], lat[:,:-1]).mean()
        dy = latlon2distance(lon[1:,:], lat[1:,:], lon[:-1,:], lat[:-1,:]).mean()
        dxdy = (dx, dy)

    if dr is None:
        dr = max(dxdy)

    # find appropriate `maxdist`
    n = 7
    slc = (slice(n, -n), slice(n, -n))
    dist_x = latlon2distance(clon, lat[slc], lon[slc], lat[slc])
    dist_y = latlon2distance(lon[slc], clat, lon[slc], lat[slc])
    _maxdist_tmp = min(np.min(dist_x[:,[0,-1]]), np.min(dist_y[[0,-1],:]))

    if maxdist is None:
        maxdist = _maxdist_tmp
    elif maxdist > _maxdist_tmp:
        warnings.warn(f"`maxdist` may be too large. Replace with {_maxdist_tmp:.2f}.")
        maxdist = _maxdist_tmp

    # find `mindist`
    if mindist is None:
        mindist = 2 * dr
    if mindist >= maxdist:
        raise ValueError("`mindist` >= `maxdist`.")
        
    # calculate axisymmetric wind speed
    radius = np.arange(mindist, maxdist, dr)
    ws_sym = CircularAvg(lon, lat, clon, clat, radius, dxdy=dxdy, intp=intp)(ws)
    
    if all(np.isnan(ws_sym)):
        msg = "Error occurs during interpolation procedure. Interpolated axisymmetric wind speed is all of NaN."
        raise ValueError(msg)
    
    return radius[np.nanargmax(ws_sym)]


class Axisymmetricity:
    """
    Compute axisymmetricity based on Miyamoto and Takemi (2013).
    See Miyamoto and Takemi (2013) for the definition.
    
    Reference
    ---------
    [1] Yoshiaki Miyamoto and Tetsuya Takemi: "A Transition Mechanism for the Spontaneous 
        Axisymmetric Intensification of Tropical Cyclones"
        J. Atmos. Sci, 70, 112-129
        https://doi.org/10.1175/JAS-D-11-0285.1
    """
    
    def __init__(self, lon, lat, clon, clat, radius, integ='trapz', dxdy=None, intp=None):
        """
        Initialization.

        Parameter
        ---------
        lon, lat : 2d array, shape = (ny, nx)
            Longtitude / latitude
        clon, clat: scaler
            The center longtitude / latitude coordinate of TC.
        radius : 1d array, shape = (nradius,)
            The radius used in calculating axisymmetricity. Unit is km.
        dxdy: Tuple(scalar, scalar). Optional
            Spatial resolution, (dx, dy). Unit is km.
            Default is None, and it would be automatically derived based on `lon` and `lat`.
        integ: str, {'trapz', 'simps'}. Optional
            Numerical integration method. 'trapz' is trapezoidal method, and `simps` is Simpson’s
            method. Default is 'trapz'.
            See scipy documents for these two methods: https://reurl.cc/X6KpYD
        intp: str, 'griddata' or 'fortran'. Optional
            Choose the interpolation method, the first one is based on `hurricane_tools.interpolate.FastGriddata`
            and the second one is based on the fortran subroutine.
            See `hurricane_tools.coord_transform.XY2RT` for more details.
        """
        if integ == 'trapz':
            self.int_func = trapz
        elif integ == 'simps':
            self.int_func = simps
        else:
            raise ValueError(f'Unavailable `integ`: {integ}. It should be "trapz" or "simps".')
            
        if dxdy is None:
            dx = latlon2distance(lon[0,0], lat[0,0], lon[0,1], lat[0,0])
            dy = latlon2distance(lon[0,0], lat[0,0], lon[0,0], lat[1,0])
            dxdy = (dx, dy)
            
        theta = np.deg2rad(np.arange(360))
        self.dtheta = theta[1] - theta[0]
        self._interp_obj = XY2RT(lon, lat, clon, clat, radius, theta, dxdy, intp)
    
    def __call__(self, var):
        """
        Calculate axisymmetricity
        
        Parameter
        ---------
        var: 2d array, shape = (ny, nx)
            The variable used to calculate axisymmetricity.
            Its shape must be equal to `lon` and `lat`.
            
        Return
        ------
        1d array, shape = (nradius,). 
            The axisymmetricity at the given radius.
        """
        dtheta = self.dtheta
        int_func = self.int_func
        
        var_rt = self._interp_obj(var)   # (nradius, ntheta)
        var_sym = var_rt.mean(axis=1)    # (nradius,)
        var_dev = var_rt - var_sym[:,np.newaxis]   # (nradius, ntheta)
        res = var_sym**2 / (var_sym**2 + int_func(var_dev**2, dx=dtheta, axis=1) / (2*np.pi))
        return res
    
    
class Rotate:
    """
    Rotate clockwise.
    """
    
    def __init__(self, lon, lat, clon, clat, rot_theta, maxdist=None, dxdy=None, intp=None):
        """
        Initialization
        
        Parameter
        ---------
        lon, lat: 2-d array, shape = (ny, nx)
            The longitude/latitude coordinates
        clon, clat: scalar
            Center coordinate
        rot_theta: scalar
            The rotate angle. Unit is radian.
            Its value must be greater than 0.
        maxdist: scalar. Optional
            The maximum search distance. Unit is km.
        mindist: scalar. Optional
            The minimum search distance. Unit is km.
        dxdy: Tuple(scalar, scalar). Optional
            Spatial resolution, (dx, dy). 
            Default is None, and it would be automatically derived based on `lon` and `lat`.
        intp: str, 'griddata' or 'fortran'. Optional
            Choose the interpolation method, the first one is based on `hurricane_tools.interpolate.FastGriddata`
            and the second one is based on the fortran subroutine.
            See `hurricane_tools.coord_transform.XY2RT` for more details.
        """
        _maxdist_tmp = self._est_maxdist(lon, lat, clon, clat)
        
        if maxdist is None:
            maxdist = _maxdist_tmp
        elif maxdist > _maxdist_tmp:
            warnings.warn(f"`maxdist` may be too large. Replace with {_maxdist_tmp:.2f}.")
            maxdist = _maxdist_tmp
            
        if dxdy is None:
            dx = latlon2distance(lon[:,1:], lat[:,1:], lon[:,:-1], lat[:,:-1]).mean()
            dy = latlon2distance(lon[1:,:], lat[1:,:], lon[:-1,:], lat[:-1,:]).mean()
            dr = max(dx, dy)
        else:
            dr = max(dxdy)
            
        if rot_theta < 0:
            raise ValueError("`rot_theta` must be greater than 0.")
        else:
            rot_theta_deg = np.rad2deg(rot_theta)
            
        radius = np.arange(0, maxdist, dr)
        theta = np.deg2rad(np.arange(360))
        theta_rot = np.deg2rad(np.arange(rot_theta_deg, rot_theta_deg+360))
        
        self._xy2rt_obj = XY2RT(lon, lat, clon, clat, radius, theta, dxdy=dxdy, intp=intp)
        self._rt2xy_obj = RT2XY(lon, lat, clon, clat, radius, theta_rot, intp=intp)
        
    def __call__(self, var):
        """
        Rotate `var` clockwise.
        
        Paramter
        --------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Rotated result. shape = (ny, nx)
        """
        return self._rt2xy_obj(self._xy2rt_obj(var))
    
    def _est_maxdist(self, lon, lat, clon, clat):
        n = 7
        slc = (slice(n, -n), slice(n, -n))
        dist_x = latlon2distance(clon, lat[slc], lon[slc], lat[slc])
        dist_y = latlon2distance(lon[slc], clat, lon[slc], lat[slc])
        _maxdist_tmp = min(np.min(dist_x[:,[0,-1]]), np.min(dist_y[[0,-1],:]))
        return _maxdist_tmp
    
    