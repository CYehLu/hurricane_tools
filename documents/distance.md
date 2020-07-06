# distance.py
<span style="color:#a77864">**latlon2distance**</span>**(lon1, lat1, lon2, lat2)**

calculate the distance (km) of two positions


******
<span style="color:#a77864">**find_lonlat_with_distance**</span>**(clon, clat, distance, xy)**

    Find the longtitude/latitude which the horizontal or vertical distance
    to the (clon, clat) is equal to 'distance'.
    
    Parameters:
    ----------
    clon, clat: scaler. The reference coordinate.
    distance: scaler. The distance (km) between target coordinate and (clon, clat).
    xy: str, 'x' or 'y'. 
        If x, it will find the coordinate which distance along latitude line is
        equal to 'distance' (The horizontal distance).
        If y, it will find the coordinate which distance along longtitude line is
        equal to 'distance' (The vertical distance).
        
    Return:
    ------
    scaler, the difference of lon/lat degree



******