from types import MethodType

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


__all__ = [
    'CWBcmapDBZ',
    'CWBcmapRain'
]


class _HybridMethod:
    """
    Make method can be called by instance or class.
    Used in `CWBcmapDBZ.show_cmap()` and `CWBcmapRain.show_cmap()`
    
    Reference : https://reurl.cc/d5l6Ag
    juanpa.arrivillaga's answer
    """
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, cls=None):
        if obj is None:
            return MethodType(self.f, cls)
        else:
            return MethodType(self.f, obj)
        

class CWBcmapDBZ:
    """
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
    """
    
    _colors_blue = [
        '#00FFFF',    # for dbz = 0 ~ 1
        '#00ECFF',    # 1 ~ 2
        '#00DAFF',
        '#00C8FF',
        '#00B6FF',
        '#00A3FF',
        '#0091FF',
        '#007FFF',
        '#006DFF',
        '#005BFF',
        '#0048FF',
        '#0036FF',
        '#0024FF',
        '#0012FF',
        '#0000FF'     # 14 ~ 15
    ]
    
    _colors_green = [
        '#00FF00',    # 15 ~ 16
        '#00F400',
        '#00E900',
        '#00DE00',
        '#00D300',
        '#00C800',
        '#00BE00',
        '#00B400',
        '#00AA00',
        '#00A000',
        '#009600',
        '#33AB00',
        '#66C000',
        '#99D500',
        '#CCEA00'    # 29 ~ 30
    ]
    
    _colors_yellow = [
        '#FFFF00',    # 30 ~ 31
        '#FFF300',
        '#FFE900',
        '#FFDE00',
        '#FFD300',
        '#FFC800',
        '#FFB800',
        '#FFA800',
        '#FF9800',
        '#FF8800',
        '#FF7800',
        '#FF6000'    # 41 ~ 42
    ]
    
    _colors_red = [
        '#FF4800',   # 42 ~ 43
        '#FF3000',
        '#FF1800',
        '#FE0000',
        '#F40000',
        '#E90000',
        '#DE0000',
        '#D30000',
        '#C80000',
        '#BE0000',
        '#B40000',
        '#AA0000',
        '#A00000',
        '#960000',
        '#AB0033'    # 56 ~ 57
    ]

    _colors_purple = [
        '#C00066',    # 57 ~ 58
        '#D50099',
        '#EA00CC',
        '#FF00FF',
        '#EA00FF',
        '#D500FF',
        '#C000FF',
        '#AB00FF',
        '#9600FF'    # 65 ~ 66
    ]
    
    colors = _colors_blue + _colors_green + _colors_yellow + _colors_red + _colors_purple
    cmap = matplotlib.colors.ListedColormap(colors)
    levels = np.arange(0, 66+1, 1)
    norm = matplotlib.colors.BoundaryNorm(levels, len(colors))
    kwargs = {'norm': norm, 'levels': levels, 'cmap': cmap}
    
    def __init__(self, n=None):
        """
        Choose coaser levels.
        
        Parameter
        ---------
        n : int, n >= 1.
            Use coaser `n` times levels.
            Default is 1 (identical to original one).
        """
        if n is None:
            n = 1
                    
        self.colors = (
            self._colors_blue[::n] 
            + self._colors_green[::n] 
            + self._colors_yellow[::n] 
            + self._colors_red[::n] 
            + self._colors_purple[::n]
        )
        self.cmap = matplotlib.colors.ListedColormap(self.colors)
        self.levels = np.arange(0, 66+1, n)
        self.norm = matplotlib.colors.BoundaryNorm(self.levels, len(self.colors))
        self.kwargs = {'norm': self.norm, 'levels': self.levels, 'cmap': self.cmap}
        
    @_HybridMethod
    def show_cmap(self, tickintv=1):
        """
        Show current color map.
        
        The xticks of plotting are `self.levels`. When `self.levels` contains non-integer,
        only one decimal place will be displayed.
        
        Parameter
        ---------
        tickitnv : int, optional
            The xtick interval of plotting. Default is 1.
        """       
        levels = self.levels
        cmap = self.cmap
        
        plt.figure(figsize=(8, 24))
        plt.imshow([np.arange(len(levels)-1)], cmap=cmap)

        xticks = np.arange(-0.5, len(levels)-0.5, 1)
        xticklabels = [f'{int(lev):d}' if int(lev) == lev else f'{lev:.1f}' for lev in levels]
        plt.xticks(xticks[::tickintv], xticklabels[::tickintv])
        plt.yticks([], [])
        plt.show()


class CWBcmapRain:
    """
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
    """
    
    colors = [
        '#CACACA',
        '#9DFEFF',
        '#01D1FD',
        '#00A6FD',
        '#0177FD',
        '#279F1A',
        '#01FB30',
        '#FFFE32',
        '#FFD328',
        '#FFA71F',
        '#FF2B06',
        '#DA2304',
        '#AB1903',
        '#AB21A3',
        '#DC2DD2',
        '#FB39FA',
        '#FFD6FE'
    ]
    cmap = matplotlib.colors.ListedColormap(colors)
    levels = np.array([0.5, 1, 2, 6, 10, 15, 20, 30, 40, 50, 70, 90, 110, 130, 150, 200, 300, 400])
    norm = matplotlib.colors.BoundaryNorm(levels, len(colors))
    kwargs = {'norm': norm, 'levels': levels, 'cmap': cmap}
    
    def __init__(self, interval=None, timeunit=None):
        """
        Choose small or large interval, daily or hourly levels.
        
        Parameter
        ---------
        interval : {'small', 'large'}
            Using small interval or large interval levels.
            Default is 'small'.
        timeunit : {'daily', 'hourly'}
            Using daily or hourly accumulated rainfall levels.
            Default is 'daily'.
        """
        if interval is None:
            interval = 'small'
        if timeunit is None:
            timeunit = 'daily'
        
        if interval == 'small' and timeunit == 'daily':
            levels = np.array([0.5, 1, 2, 6, 10, 15, 20, 30, 40, 50, 70, 90, 110, 130, 150, 200, 300, 400])
        elif interval == 'large' and timeunit == 'daily':
            levels = np.array([0.5, 10, 20, 60, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 1900])
        elif interval == 'small' and timeunit == 'hourly':
            levels = np.array([0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17])
        elif interval == 'large' and timeunit == 'hourly':
            levels = np.array([0.5, 1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])
            
        self.levels = levels
        self.norm = matplotlib.colors.BoundaryNorm(self.levels, len(self.colors))
        self.kwargs = {'norm': self.norm, 'levels': self.levels, 'cmap': self.cmap}
        
    def set_levels(self, levels):
        """
        Set new levels. `norm` and `kwargs` will also change.
        
        Note : levels.size should equal to 18, because len(self.colors) = 17.
        """
        self.levels = levels
        self.norm = matplotlib.colors.BoundaryNorm(self.levels, len(self.colors))
        self.kwargs = {'norm': self.norm, 'levels': self.levels, 'cmap': self.cmap}
        
    @_HybridMethod
    def show_cmap(self, tickintv=1):
        """
        Show current color map.
        
        The xticks of plotting are `self.levels`. When `self.levels` contains non-integer,
        only one decimal place will be displayed.
        
        Parameter
        ---------
        tickitnv : int, optional
            The xtick interval of plotting. Default is 1.
        """       
        levels = self.levels
        cmap = self.cmap
        
        plt.figure(figsize=(8, 24))
        plt.imshow([np.arange(len(levels)-1)], cmap=cmap)

        xticks = np.arange(-0.5, len(levels)-0.5, 1)
        xticklabels = [f'{int(lev):d}' if int(lev) == lev else f'{lev:.1f}' for lev in levels]
        plt.xticks(xticks[::tickintv], xticklabels[::tickintv])
        plt.yticks([], [])
        plt.show()
