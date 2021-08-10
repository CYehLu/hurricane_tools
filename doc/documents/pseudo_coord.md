# pseudo_coord  

[[source](../.././hurricane_tools//pseudo_coord.py)]  

<span style="color:#a77864">**lonlat2xy**</span>**(lon, lat, clon, clat)**

    Convert to distance-based coordinate.
    
    Parameter
    ---------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and latitude
    clon, clat : scalar
        Center longtitude-latitude coordinate
        
    Return
    ------
    X, Y : 2d array, shape = (ny, nx)
        `X` is the zonal distance between each lon-lat coordinate points and
        the center coordinate, and `Y` is the meridional distance. 
        Unit is km.



******
<span style="color:#a77864">**xy2lonlat**</span>**(X, Y, clon=None, clat=None, equal_dist=True)**

    Convert to longtitude-latitude coordinate.
    
    Parameter
    ---------
    X, Y : 2d array, shape = (ny, nx)
        `X` is the zonal distance between each grid points and the center
        coordinate (`clon`, `clat`), and `Y` is the meridional distance.
    clon, clat : scalar, optional
        Center longtitude-latitude coordinate.
        Default is clon=120, clat=20
    equal_dist : bool, optional
        Determine if the zonal/meridional distance between each zonal/meridional
        grids is nearly equal.
        This can only be True now.
        
    Return
    ------
    lon, lat : 2d array, shape = (ny, nx)
        Longtitude and Latitude.



******