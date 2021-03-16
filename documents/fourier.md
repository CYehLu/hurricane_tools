# fourier  

[[source](.././hurricane_tools//fourier.py)]  

<span style="color:#a77864">**interp_xy_closure**</span>**(radius, theta, x, y, center=None)**

    Same as `interp_xy`, but return a closure function to perform the interpolation.



******
<span style="color:#a77864">**interp_xy**</span>**(radius, theta, values, x, y, center=None)**

    Interpolating data from radius-theta coordinate to x-y coordinate.



******
class <span style="color:#a77864">**CircularFourier**</span>

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



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **circular\_fft** </font> | Perform FFT on circles. |
| <font color="#a77864"> **decomp** </font> | Decompute variable and extract the wave number `ns` components respectively. |


<span style="color:#cca99b">CircularFourier</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, dxdy, dtheta=None)**

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

  
<span style="color:#cca99b">CircularFourier</span>.<span style="color:#a77864">**circular\_fft**</span>**(self, var)**

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