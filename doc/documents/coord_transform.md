# coord_transform  

[[source](../.././hurricane_tools//coord_transform.py)]  

class <span style="color:#a77864">**XY2RT**</span>

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



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **\_\_call\_\_** </font> | Interpolate data from cartesian (x-y) to polar (radius-theta) coordinate. |


<span style="color:#cca99b">XY2RT</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, theta=None, dxdy=None, intp=None)**

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

  
<span style="color:#cca99b">XY2RT</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, var)**

        Interpolate data from cartesian (x-y) to polar (radius-theta) coordinate.
        
        Parameter
        ---------
        var: 2-d array, shape = (ny, nx)
            Its shape should equal to `lon` and `lat`.
        
        Return
        ------
        Interpolating result, shape = (nradius, ntheta)

  
******
class <span style="color:#a77864">**RT2XY**</span>

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



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. |
| <font color="#a77864"> **\_\_call\_\_** </font> | Interpolate data from polar (radius-theta) to cartesian (x-y) coordinate. |


<span style="color:#cca99b">RT2XY</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, lon, lat, clon, clat, radius, theta, intp=None)**

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

  
<span style="color:#cca99b">RT2XY</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, var)**

        Interpolate data from polar (radius-theta) to cartesian (x-y) coordinate.
        
        Parameter
        ---------
        var: 2-d array, shape = (nradius, ntheta)
        
        Return
        ------
        Interpolating result, shape = (ny, nx)

  
******
class <span style="color:#a77864">**Interpz3d**</span>

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



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialize with pressure and levels. |
| <font color="#a77864"> **interp** </font> | Interpolate variable from the original vertical coordinate to `zdata` coordinate. |


<span style="color:#cca99b">Interpz3d</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, zdata, level, missing_value=np.nan)**

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

  
<span style="color:#cca99b">Interpz3d</span>.<span style="color:#a77864">**interp**</span>**(self, \*var)**

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

  
******