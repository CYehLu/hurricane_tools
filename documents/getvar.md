# getvar  

[[source](.././hurricane_tools//getvar.py)]  

class <span style="color:#a77864">**GetVar**</span>

    Get variables. It is similar to `wrf.getvar` by wrf-python,
    but here I store every intermediate variables to speed up.



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization the instance. |
| <font color="#a77864"> **close** </font> | Close attribute of GetVar instance. |
| <font color="#a77864"> **get** </font> | Get variable by its name. |


<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, filename, timeidx=None)**

        Initialization the instance.

        Parameters
        ----------
        filename : str
            The netCDF file name
        timeidx : optional, int or slice()
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

  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**get**</span>**(self, var_name)**

        Get variable by its name.

        Parameters
        ----------
        var_name : str
            Variable name.
            It can be the variable name of the netCDF file, or some diagnosis variables
            list below:
                'slp'  --  Sea Level Pressure
                'tk'   --  Temperature (unit: K)

  
******