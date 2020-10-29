# plot  

[[source](.././hurricane_tools//plot.py)]  

class <span style="color:#a77864">**CWBcmapDBZ**</span>

    Colormap, colors, contour levels and norm of CWB dbz pictures.
    
    Example
    -------
    >>> dbz = _generate_fake_data()
    >>> # 1. default setting
    >>> plt.contourf(
    ...     dbz, 
    ...     norm=CWBcmapDBZ.norm, 
    ...     levels=CWBcmapDBZ.levels, 
    ...     cmap=CWBcmapDBZ.cmap
    ... )
    >>> plt.show()
    >>> 
    >>> # 2. more coarse color levels
    >>> cwbdbz = CWBcmapDBZ(n=4)    # 4 times coarser
    >>> plt.contourf(
    ...     dbz, 
    ...     norm=cwbdbz.norm, 
    ...     levels=cwbdbz.levels, 
    ...     cmap=cwbdbz.cmap
    ... )
    >>> plt.show()
    
    Attributes 
    ----------
    colors : List[str]
        List of color hex. default length = 66
    cmap : matplotlib.colors.ListedColormap
        Colormap
    levels : np.array, default is [0, 1, ..., 66].
        Contour levels
    norm : matplotlib.colors.BoundaryNorm
        Can be used in `plt.contourf`



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Choose coaser levels. |


<span style="color:#cca99b">CWBcmapDBZ</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, n=None)**

        Choose coaser levels.
        
        Parameter
        ---------
        n : int, n >= 1.
            Use coaser `n` times levels.
            Default is 1 (identical to original one).

  
******
class <span style="color:#a77864">**CWBcmapRain**</span>

    Colormap, colors, contour levels and norm of CWB accumulated daily/hourly rainfall pictures.
    
    There are 4 options : 
        interval = 'small' or 'large', and time unit = 'daily' or 'hourly'
    They only differ in `levels` and `norm`. The max of levels for 4 options :
        | ----------------------- |
        |         daily  hourly   |
        | small     400      17   |
        | large    1900      75   |
        | ----------------------- |
    Default is 'small' and 'daily'.
    
    Example
    -------
    >>> rain = _generate_fake_data()
    >>> # 1. default setting
    >>> plt.contourf(
    ...     rain, 
    ...     norm=CWBcmapRain.norm, 
    ...     levels=CWBcmapRain.levels, 
    ...     cmap=CWBcmapRain.cmap
    ... )
    >>> plt.show()
    >>> 
    >>> # 2. choose large interval and hourly time unit
    >>> cwbrain = CWBcmapRain(interval='large', timeunit='hourly')
    >>> plt.contourf(
    ...     rain, 
    ...     norm=cwbrain.norm, 
    ...     levels=cwbrain.levels, 
    ...     cmap=cwbrain.cmap
    ... )
    >>> plt.show()
    
    Attributes 
    ----------
    colors : List[str]
        List of color hex. length = 17
    cmap : matplotlib.colors.ListedColormap
        Colormap
    levels : np.array
        Contour levels
        For small interval (daily), max of levels is 400. And 1900 for large interval.
    norm : matplotlib.colors.BoundaryNorm
        Can be used in `plt.contourf`



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Choose small or large interval, daily or hourly levels. |


<span style="color:#cca99b">CWBcmapRain</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, interval=None, timeunit=None)**

        Choose small or large interval, daily or hourly levels.
        
        Parameter
        ---------
        interval : {'small', 'large'}
            Using small interval or large interval levels.
            Default is 'small'.
        timeunit : {'daily', 'hourly'}
            Using daily or hourly accumulated rainfall levels.
            Default is 'daily'.

  
******