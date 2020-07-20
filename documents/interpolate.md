# interpolate
class <span style="color:#a77864">**FastGriddata**</span>

    Faster version of scipy.interpolate.griddata (in repeatly interpoting 
    case).
    
    When performs interpolation repeatly, and parameters remain the same
    except the interpolated values, this can reduce some computational 
    time by reuse the "triangle grid" during the process of griddata.
    
    Codes are modified from : https://reurl.cc/xZ50Yb
    The first and the second answers.
    
    Example:
    -------
    >>> import time
    >>> from scipy.interpolate import griddata
    >>> from hurricane_tools.interpolate import FastGriddata
    >>>
    >>> # construct data points
    ... X_, Y_ = np.meshgrid(np.linspace(-10, 10), np.linspace(-10, 10))
    ... xy = np.vstack((X_.ravel(), Y_.ravel())).T   # shape=(100*100, 2)
    ... uv = np.random.rand(1000, 2) * 20 - 10
    >>>
    >>> # construct 3 different interpolating values
    ... vals1 = np.ravel(np.cos(X_) * np.cos(Y_))
    ... vals2 = np.ravel(np.sin(X_) + 2 * np.cos(Y_))
    ... vals3 = np.ravel(np.cos(X_**2) * np.cos(Y_))
    >>>
    >>> # using scipy.interpolate.gridata -- it is slower
    ... time_start = time.time()
    ... res1_g = griddata(xy, vals1, uv)
    ... res2_g = griddata(xy, vals2, uv)
    ... res3_g = griddata(xy, vals3, uv)
    ... time_end = time.time()
    ... print('by using griddata :', time_end - time_start)   # 0.1519 sec
    by using griddata : 0.15190553665161133
    >>>
    >>> # using FasterGriddata -- faster than apply griddata 3 times
    ... time_start = time.time()
    ... interp_obj = FastGriddata(xy, uv)
    ... res1_f = interp_obj.interpolate(vals1)
    ... res2_f = interp_obj.interpolate(vals2)
    ... res3_f = interp_obj.interpolate(vals3)
    ... time_end = time.time()
    ... print('by using FastGriddata :', time_end - time_start)   # 0.0543 sec
    by using FastGriddata : 0.05439400672912598
    >>> 
    >>> # check results
    ... np.allclose(res1_g, res1_f)
    True
    >>> np.allclose(res2_g, res2_f)
    True
    >>> np.allclose(res3_g, res3_f)
    True

| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **__init__** </font> | xy: grid points, shape = (n, 2) uv: interpolated points, shape = (m, 2) |
| <font color="#a77864"> **__call__** </font> |  |
| <font color="#a77864"> **interpolate** </font> | values: values on `xy` coordinates. shape = (n, 2) return interpolating result, shape = (m, 2) |


<span style="color:#cca99b">FastGriddata</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, xy, uv, d=2)**

        xy: grid points, shape = (n, 2)
        uv: interpolated points, shape = (m, 2)

  
<span style="color:#cca99b">FastGriddata</span>.<span style="color:#a77864">**\_\_call\_\_**</span>**(self, values)**


  
<span style="color:#cca99b">FastGriddata</span>.<span style="color:#a77864">**interpolate**</span>**(self, values)**

        values: values on `xy` coordinates. shape = (n, 2)
        return interpolating result, shape = (m, 2)

  
******