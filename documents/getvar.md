# getvar  

[[source](.././hurricane_tools//getvar.py)]  

class <span style="color:#a77864">**GetVar**</span>

    Get variables. It is similar to `wrf.getvar` by wrf-python,
    but here I store every intermediate variables to speed up.



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Initialization. `filename` should be a str. |
| <font color="#a77864"> **close** </font> | close netCDF file instance |
| <font color="#a77864"> **get** </font> | Get variable by its name. |


<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, filename)**

    Initialization. `filename` should be a str.
  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**close**</span>**(self)**

    close netCDF file instance
  
<span style="color:#cca99b">GetVar</span>.<span style="color:#a77864">**get**</span>**(self, var_name)**

    Get variable by its name.
  
******