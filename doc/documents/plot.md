# plot  

[[source](../.././hurricane_tools//plot.py)]  

class <span style="color:#a77864">**CWBcmapDBZ**</span>

    Colormap, colors, contour levels and norm of CWB dbz pictures.
    
    Example
    -------
    >>> dbz = _generate_fake_data()
    >>> # 1. default setting
    >>> plt.contourf(dbz, **CWBcmapDBZ.kwargs)
    >>> plt.show()
    >>> 
    >>> # 2. coarser color levels
    >>> cwbdbz = CWBcmapDBZ(n=4)    # 4 times coarser
    >>> plt.contourf(dbz, **cwbdbz.kwargs)
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
    kwargs : Dict
        kwargs = {'norm': ..., 'cmap': ..., 'levels': ...}
        Using `**CWBcmapDBZ.kwargs` in `contourf`
            >>> plt.contourf(dbz, **CWBcmapDBZ.kwargs)
        is equivalent to 
            >>> plt.contourf(
            ...     dbz, 
            ...     levels=CWBcmapDBZ.levels
            ...     cmap=CWBcmapDBZ.cmap,
            ...     norm=CWBcmapDBZ.norm
            ... )



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Choose coaser levels. |
| <font color="#a77864"> **show\_colorbar** </font> | Show current colorbar. |


<span style="color:#cca99b">CWBcmapDBZ</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, n=None)**

        Choose coaser levels.
        
        Parameter
        ---------
        n : int, n >= 1.
            Use coaser `n` times levels.
            Default is 1 (identical to original one).

  
<span style="color:#cca99b">CWBcmapDBZ</span>.<span style="color:#a77864">**show\_colorbar**</span>**(self, tickintv=None, ticks=None, figsize=None)**

        Show current colorbar.
        
        The ticks of colorbar are `self.levels`. When `self.levels` contains non-integer,
        only one decimal place will be displayed.
        
        Parameter
        ---------
        tickintv : int, optional
            The ticks interval of plotting. Default is 1.
        ticks : array-like, optional
            The ticks of colorbar. 
            `tickintv` has priority over `ticks` when both of them are given.
        figsize : Tuple(int, int). optional
            Figure size. Default is (10, 1).
            
        Return
        ------
        Instance of `matplotlib.colorbar.ColorbarBase`

  
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
    >>> plt.contourf(rain, **CWBcmapRain.kwargs)
    >>> plt.show()
    >>> 
    >>> # 2. choose large interval and hourly time unit
    >>> cwbrain = CWBcmapRain(interval='large', timeunit='hourly')
    >>> plt.contourf(rain, **cwbrain.kwargs)
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
    kwargs : Dict
        kwargs = {'norm': ..., 'cmap': ..., 'levels': ...}
        Using `**CWBcmapRain.kwargs` in `contourf`
            >>> plt.contourf(rain, **CWBcmapRain.kwargs)
        is equivalent to 
            >>> plt.contourf(
            ...     rain, 
            ...     levels=CWBcmapRain.levels
            ...     cmap=CWBcmapRain.cmap,
            ...     norm=CWBcmapRain.norm
            ... )



| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Choose small or large interval, daily or hourly levels. |
| <font color="#a77864"> **set\_levels** </font> | Set new levels. `norm` and `kwargs` will also change. |
| <font color="#a77864"> **show\_colorbar** </font> | Show current colorbar. |


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

  
<span style="color:#cca99b">CWBcmapRain</span>.<span style="color:#a77864">**set\_levels**</span>**(self, levels)**

        Set new levels. `norm` and `kwargs` will also change.
        
        Note : levels.size should equal to 18, because len(self.colors) = 17.

  
<span style="color:#cca99b">CWBcmapRain</span>.<span style="color:#a77864">**show\_colorbar**</span>**(self, tickintv=None, ticks=None, figsize=None)**

        Show current colorbar.
        
        The ticks of colorbar are `self.levels`. When `self.levels` contains non-integer,
        only one decimal place will be displayed.
        
        Parameter
        ---------
        tickintv : int, optional
            The ticks interval of plotting. Default is 1.
        ticks : array-like, optional
            The ticks of colorbar. 
            `tickintv` has priority over `ticks` when both of them are given.
        figsize : Tuple(int, int). optional
            Figure size. Default is (10, 1).
            
        Return
        ------
        Instance of `matplotlib.colorbar.ColorbarBase`

  
******