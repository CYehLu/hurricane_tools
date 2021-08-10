# distance  

[[source](../.././hurricane_tools//distance.py)]  

<span style="color:#a77864">**latlon2distance**</span>**(lon1, lat1, lon2, lat2)**

    Calculate the distance (km) between two points in longitude/latitude coordinate.
    Units of parameters are degree.



******
<span style="color:#a77864">**find_dlonlat_by_distance**</span>**(clon, clat, distance, xy)**

    Find the difference of longitude/latitude which the zonal/meridional
    distance to the (clon, clat) is equal to `distance`.
    
    Parameters
    ----------
    clon, clat: scaler, unit: degree
        The reference coordinate 
    distance: scaler, unit: km
        The distance between target coordinate and (clon, clat).
    xy: str, 'x' or 'y'. 
        If xy == `x`:
            return value = `dlon`, that is the distance between (clon, clat) and
            (clon+dlon, clat) = `distance`.
        If xy == `y`:
            return value = `dlat`, that is the distance between (clon, clat) and
            (clon, clat+dlat) = `distance`.
        
    Return
    ------
    scaler (unit: degree), the difference of longitude or latitude degree.



******