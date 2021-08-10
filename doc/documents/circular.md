# circular  

[[source](../.././hurricane_tools//circular.py)]  

class <span style="color:#a77864">**CircularAvg**</span>

    Compute the azithmual mean field (axisymmetric field).



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **\_\_call\_\_** </font> | Calculate the azimuthal mean. |


<span style="color:#cca99b">CircularAvg</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None)**

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

  
<span style="color:#cca99b">CircularAvg</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, values)**

        Calculate the azimuthal mean.
        
        Parameter
        ---------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Azimuthal mean, shape = (nradius,)

  
******
class <span style="color:#a77864">**CircularFourier**</span>

    Fourier transform in cylindrical coordinate.
    
    Example
    -------
    >>> lon, lat, clon, clat, qvapor = ...   # get data
    >>> lon.shape     # (ny, nx)
    (201, 201)
    >>> radius = np.arange(0, 100)   # unit: km
    >>> cfourier = CircularFourier(lon, lat, clon, clat, radius)
    >>> k, qv_hat = cfourier.circular_fft(qvapor)   # spectral of `qvapor` on difference radius
    >>>
    >>> qv_wns = cfourier.decomp(qvapor, [0, 1, 2, 3])   # decompose `qvapor` to wavenumber-0 to 3 components
    >>> len(qv_wns)       # List[Array]
    4
    >>> qv_wns[0].shape   # (ny, nx)
    (201, 201)



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **circular\_fft** </font> | Perform FFT in cylindrical coordinate. |
| <font color="#a77864"> **decomp** </font> | Decompute variable and extract the wave number `ns` components respectively. |


<span style="color:#cca99b">CircularFourier</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, dxdy=None, intp=None)**

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

  
<span style="color:#cca99b">CircularFourier</span>.<span style="color:#a77864">**circular\_fft**</span>**(self, var)**

        Perform FFT in cylindrical coordinate.
        
        Parameter
        ---------
        var : 2-d array, shape = (ny, nx)
            Variable to perform FFT
            
        Return
        ------
        k : 1-d array, shape = (m,)
            Wave numbers
        var_hat : 2-d array, shape = (nradius, m)
            FFT spectral density

  
<span style="color:#cca99b">CircularFourier</span>.<span style="color:#a77864">**decomp**</span>**(self, var, ns)**

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

  
******
<span style="color:#a77864">**rmw**</span>**(lon, lat, ws, clon, clat, maxdist=None, mindist=None, dr=None, dxdy=None, intp=None)**

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



******
class <span style="color:#a77864">**Axisymmetricity**</span>

    Compute axisymmetricity based on Miyamoto and Takemi (2013).
    See Miyamoto and Takemi (2013) for the definition.
    
    Reference
    ---------
    [1] Yoshiaki Miyamoto and Tetsuya Takemi: "A Transition Mechanism for the Spontaneous 
        Axisymmetric Intensification of Tropical Cyclones"
        J. Atmos. Sci, 70, 112-129
        https://doi.org/10.1175/JAS-D-11-0285.1



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **\_\_call\_\_** </font> | Calculate axisymmetricity |


<span style="color:#cca99b">Axisymmetricity</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, integ='trapz', dxdy=None, intp=None)**

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

  
<span style="color:#cca99b">Axisymmetricity</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, var)**

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

  
******
class <span style="color:#a77864">**Rotate**</span>

    Rotate clockwise.



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization |
| <font color="#a77864"> **\_\_call\_\_** </font> | Rotate `var` clockwise. |


<span style="color:#cca99b">Rotate</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, rot_theta, maxdist=None, dxdy=None, intp=None)**

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

  
<span style="color:#cca99b">Rotate</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, var)**

        Rotate `var` clockwise.
        
        Paramter
        --------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Rotated result. shape = (ny, nx)

  
******