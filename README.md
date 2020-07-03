center.py
------
| Function | Description |
| -------- | ----------- |
| <font color="#a77864"> **ty_center_mslp** </font> | Find typhoon center by minimum sea level pressure grid (no weighted). |
| <font color="#a77864"> **weighted_ty_center** </font> | Giving first guess typhoon center longtitude and latitude, calculate the new center by weighted mothod. |


******
parse_wrf_center_mesg.py
------
| Function | Description |
| -------- | ----------- |
| <font color="#a77864"> **parse_wrf_rsl_error** </font> | Parse `rsl.error.0000` outputed from WRF to get the center information. If parse file successfully, it would create `center.txt` in the current folder. |
| <font color="#a77864"> **find_centers_nearest_time** </font> | Find the hurricane center data at specific times from the file which contains center information based on WRF rsl.error.0000. |


******
circular.py
------
| Function | Description |
| -------- | ----------- |
| <font color="#a77864"> **circular_avg** </font> | Calculate circular mean. |
| <font color="#a77864"> **circular_avg_closure** </font> | Calculate circular mean. |
| <font color="#a77864"> **rmw** </font> | Find TC RMW |
| <font color="#a77864"> **axissymmetricity** </font> | var: shape=(ny, nx) radius: shape=(nr,) |


******
distance.py
------
| Function | Description |
| -------- | ----------- |
| <font color="#a77864"> **latlon2distance** </font> | calculate the distance (km) of two positions |
| <font color="#a77864"> **find_lonlat_with_distance** </font> | Find the longtitude/latitude which the horizontal or vertical distance to the (clon, clat) is equal to 'distance'. |


******
transform.py
------
| Function | Description |
| -------- | ----------- |
| <font color="#a77864"> **uv2vrvt_rt** </font> | u, v: (ny, nx) -> vr, vt: (nradius, ntheta) |
| <font color="#a77864"> **uv2vrvt_xy** </font> | u, v: (ny, nx) -> vr, vt: (ny, nx) |


******