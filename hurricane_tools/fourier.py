import numpy as np
from scipy.interpolate import griddata

from hurricane_tools.circular import interp_circle
from hurricane_tools.distance import latlon2distance
from hurricane_tools.interpolate import FastGriddata


__all__ = [
    'interp_xy_closure',
    'interp_xy',
    'CircularFourier'
]



def interp_xy_closure(radius, theta, x, y, center=None):
    """
    Same as `interp_xy`, but return a closure function to perform the interpolation.
    """
    if center is None:
        cx = 0
        cy = 0
    else:
        cx, cy = center
        
    if x.ndim != y.ndim:
        raise ValueError(f"x.ndim ({x.ndim}) != y.ndim ({y.ndim})")
    elif x.ndim >= 3 or y.ndim >= 3:
        raise ValueError(f"x.ndim and y.ndim should <= 2")
    else:
        if x.ndim == 1 and y.ndim == 1:
            X, Y = np.meshgrid(x, y)
        elif x.ndim == 2 and y.ndim == 2:
            X, Y = x, y
            
    t_mesh, r_mesh = np.meshgrid(theta, radius)
    polar_coord = np.vstack((r_mesh.ravel(), t_mesh.ravel())).T
    
    interp_coord_r = np.sqrt((X-cx)**2 + (Y-cy)**2)
    interp_coord_t = np.arctan2(Y-cy, X-cx)
    interp_coord_t[interp_coord_t < 0] += 2 * np.pi
    interp_coord = np.vstack((interp_coord_r.ravel(), interp_coord_t.ravel())).T
    
    fgd_obj = FastGriddata(polar_coord, interp_coord)
    
    def inner(values):
        return fgd_obj.interpolate(values.ravel()).reshape(X.shape)
    
    return inner


def interp_xy(radius, theta, values, x, y, center=None):
    """
    Interpolating data from radius-theta coordinate to x-y coordinate.
    """
    interp_func = interp_xy_closure(radius, theta, x, y, center)
    return interp_func(values)



class CircularFourier:
    """
    Perform Fourier transformation on circles.
    
    Example
    -------
    >>> # get data
    >>> import hurricane_tools as ht
    >>> gv = ht.getvar.GetVar('test.nc')
    >>> lon = gv.get('XLONG')
    >>> lat = gv.get('XLAT')
    >>> dbz = gv.get('dbz')
    >>> dx = gv.ncfile.DX / 1000
    >>> dy = gv.ncfile.DY / 1000
    >>> gv.close()
    >>>
    >>> # find TC center
    >>> clon, clat = ht.center.weighted_tc_center(lon, lat, slp)
    >>>
    >>> # Perform fourier and extract 0 ~ 3 wave number components
    >>> radius = np.arange(0, 180, 1)
    >>> cfourier = ht.fourier.CircularFourier(lon, lat, clon, clat, radius, dxdy=(dx, dy))
    >>> ns = [0, 1, 2, 3]
    >>> dbz_ns = cfourier.decomp(dbz[0,:,:], ns)    # List[Array]
    >>>
    >>> # plot
    >>> import matplotlib.pyplot as plt
    >>> plt.figure()
    >>> plt.contourf(lon, lat, dbz)    # original data
    >>> 
    >>> for n in ns:
    ...     plt.figure()
    ...     plt.contourf(lon, lat, dbz_ns[n])
    ...     plt.title(f'wave number component : {n}')
    """
    
    def __init__(self, lon, lat, clon, clat, radius, dxdy, dtheta=None):
        """
        Initialization.
        
        Parameter
        ---------
        lon, lat : 2-d array, shape = (ny, nx)
            Longtitude and Latitude
        clon, clat : scalar
            Central longtitude / latitude
        radius : 1-d array, shape = (nr,)
            The radius of circles. Unit is km.
        dxdy : 2-scalar-elements tuple
            Spatial resolution
        dtheta : scalar, optional
            Azimuthal resolution (degree). Default is 1 degree.
        """
        self.lon = lon
        self.lat = lat
        self.clon = clon
        self.clat = clat
        self.radius = radius
        self.dxdy = dxdy
        
        self._XY = None    # for `self.decomp`
        
        if dtheta is None:
            self.dtheta = np.deg2rad(1)
            self._theta = np.arange(*np.deg2rad([0, 360, 1]))
        else:
            self.dtheta = dtheta
            self._theta = np.arange(*np.deg2rad([0, 360, dtheta]))
        
    def circular_fft(self, var):
        """
        Perform FFT on circles.
        
        Parameter
        ---------
        var : 2-d array, shape = (ny, nx)
            Variable to perform FFT
            
        Return
        ------
        k : 1-d array
            Wave numbers
        var_hat : 2-d array, shape = (n_radius, n_k)
            FFT result
        """
        var_circle = interp_circle(
            self.lon, self.lat, var,
            self.clon, self.clat, self.radius, self._theta,
            dxdy=self.dxdy, coord='lonlat'
        )
        var_hat = np.fft.fft(var_circle, axis=1)    # (nr, ntheta/nk)
        k = 2 * np.pi * np.fft.fftfreq(var_circle.shape[1], d=self.dtheta)
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
        if self._XY is None:        
            lon = self.lon
            lat = self.lat
            clon = self.clon
            clat = self.clat

            # distance matrix (no directionality, all positive)
            dist_x = latlon2distance(lon, lat, clon, lat)
            dist_y = latlon2distance(lon, lat, lon, clat)

            # convert distance matrix into cartesian coordinate
            X = dist_x
            X[lon < clon] *= -1
            Y = dist_y
            Y[lat < clat] *= -1
            
            # store in instance
            self._XY = (X, Y)    
            
        else:
            X, Y = self._XY
            
        # perform fft on every radius
        k, var_hat = self.circular_fft(var)
        
        # function of convert variable from polar coord to xy coord
        polar2xy_func = interp_xy_closure(self.radius, self._theta, X, Y, center=(0, 0))
        
        res = []
        for n in ns:
            # filter out all wave number except `n`
            var_hat_tmp = np.copy(var_hat)
            var_hat_tmp[:, ~(np.abs(k) == n)] = 0
            var_n = np.fft.ifft(var_hat_tmp, axis=1).real
            var_n_xy = polar2xy_func(var_n)
            #var_n_xy = interp_xy(var_n, self.radius, self._theta, X, Y, center=(0, 0))
            res.append(var_n_xy)
            
        return res