# getvar  

[[source](../.././hurricane_tools//getvar.py)]  

<span style="color:#a77864">**destagger**</span>**(var, stagger_dim)**

    Convert variable from staggered to unstagger grid.
    
    This is almost exactly the same as `wrf.destagger`, but I simplified to 
    destagger process and restricted the `stagger_dim` to one of z, y, or, x
    dimension to improve the efficiency.
    
    Parameter
    ---------
    var : ndarray
        The variable on the staggered grid. The dimensions of this variable
        must be (..., z, y, x) or (..., y, x). 
    stagger_dim : int
        The dimension index to destagger.
        Only the rightmost 3 dimensions (z, y, or x) can be destaggered.
    
    Return
    ------
    Variable on the unstaggered grid.



******
class <span style="color:#a77864">**GetVar**</span>

    Get variables. It is similar to `wrf.getvar` by wrf-python, but here I 
    store every intermediate variables, reduce the amount of function calling
    and rewrite the fortran functions to speed up.



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization the instance. |
| <font color="#a77864"> **close** </font> | Close attribute of GetVar instance. |
| <font color="#a77864"> **get\_times** </font> | Get `Times` variable from netCDF file. |
| <font color="#a77864"> **get** </font> | Get variable by its name. |


<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, filename, timeidx=None)**

        Initialization the instance.

        Parameters
        ----------
        filename : str
            The netCDF file name
        timeidx : optional, int, list, or slice()
            The time index of variables.
            Default is None, and it would use `slice(0, None)`, the all time index.

  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**close**</span>**(self, ncfile=True, variables=False)**

        Close attribute of GetVar instance.

        Parameters
        ----------
        ncfile : optional, bool
            Close GetVar.ncfile or not. Default is True.
        variables : optional, bool
            Delete GetVar.variables or not. Default is False.

  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**get\_times**</span>**(self, dt=None)**

        Get `Times` variable from netCDF file.
        
        Directly get `Times` from netCDF file would return a (ntime, 19) array, each column
        of array is a `numpy.bytes_` and it is not easy to read.
        This method would transform this np.bytes_ array, and return a list of string (or 
        datetime object, pandas DatetimeIndex) which is more convenient to read.
        
        Paramter
        --------
        dt : str, 'python' or 'pandas'. optional
            Default is None, it would return a list of string.
            If 'python', it would return a list of datetime object.
            If 'pandas', it would return a pandas.DatatimeIndex

  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**get**</span>**(self, var_name, squeeze=True, filled=True)**

        Get variable by its name.

        Parameters
        ----------
        var_name : str
            Variable name.
            It can be the variable name of the netCDF file, or some diagnosis variables
            list below:
                'slp'   --  Sea Level Pressure
                'tk'    --  Temperature (unit: K)
                'pres'  --  Pressure (unit: hPa)
                'geopt' --  Geopotential (unit: m2 s-2)
                'theta' --  (Full) Potential Temperature (unit: K)
                'dbz'   --  Radar Reflectivity 
                'avo'   --  Absolute Vorticity (unit: 10-5 s-1)
                'pvo'   --  Potential Vorticity (unit: PVU)
                
        squeeze : bool, default is True
            Determine whether to squeeze the returned variable.
            
        filled : bool or scalar. Default is True.
            Determine the filled value of returned variable, if it is a `MaskArray`.
            
            If filled = True:
                Filled np.nan to variable. Returned variable type is `ndarray`.
            If filled = False:
                Do not fill any value to variable.
                Returned variable type is `MaskArray` (if wrf-output-variable) or `ndarray` (if 
                diagnois variable).
            If filled = scalar (int or float):
                Filled given value to variable. Returned variable type is `ndarray`.
                When filled = np.nan, it is equivalent to filled = True.

  
******