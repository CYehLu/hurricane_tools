# hurricane_tools
documents of all modules -> `./documents`


axisym_vortex
------
[document](./documents/axisym_vortex.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **rankine_vortex** </font> | Classic rankine vortex. |
| <font color="#a77864"> **holland80** </font> | The empirical formula of TC axisymmetric structure by Holland (1980). |
| <font color="#a77864"> **willoughby04** </font> | The empirical formula of TC axisymmetric structure by Willoughby (2006). |


******
bst_parser
------
[document](./documents/bst_parser.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **DownloadWarning** </font> |  |
| <font color="#a77864"> **JMAbstParser** </font> | Parse JMA best track information. |


******
center
------
[document](./documents/center.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **tc_center_mslp** </font> | Find TC center by finding the grid of minimum sea level pressure. |
| <font color="#a77864"> **weighted_tc_center** </font> | Calculate TC center by weighted method. |


******
circular
------
[document](./documents/circular.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **CircularAvg** </font> | Compute the azithmual mean field (axisymmetric field). |
| <font color="#a77864"> **CircularFourier** </font> | Fourier transform in cylindrical coordinate. |
| <font color="#a77864"> **rmw** </font> | Find the radius of maximum wind. |
| <font color="#a77864"> **Axisymmetricity** </font> | Compute axisymmetricity based on Miyamoto and Takemi (2013). See Miyamoto and Takemi (2013) for the definition. |
| <font color="#a77864"> **Rotate** </font> | Rotate clockwise. |


******
coord_transform
------
[document](./documents/coord_transform.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **XY2RT** </font> | Coordinate transformation from cartesian (x-y) to polar (radius-theta) coordinate. |
| <font color="#a77864"> **RT2XY** </font> | Coordinate transformation from polar (radius-theta) to cartesian (x-y) coordinate. |
| <font color="#a77864"> **Interpz3d** </font> | Interpolating variables on specified vertical coordinate. |


******
distance
------
[document](./documents/distance.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **latlon2distance** </font> | Calculate the distance (km) between two points in longitude/latitude coordinate. Units of parameters are degree. |
| <font color="#a77864"> **find_dlonlat_by_distance** </font> | Find the difference of longitude/latitude which the zonal/meridional distance to the (clon, clat) is equal to `distance`. |


******
getvar
------
[document](./documents/getvar.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **destagger** </font> | Convert variable from staggered to unstagger grid. |
| <font color="#a77864"> **GetVar** </font> | Get variables. It is similar to `wrf.getvar` by wrf-python, but here I store every intermediate variables, reduce the amount of function calling and rewrite the fortran functions to speed up. |


******
interpolate
------
[document](./documents/interpolate.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **FastGriddata** </font> | Faster version of scipy.interpolate.griddata (in repeatly interpolating case). |


******
plot
------
[document](./documents/plot.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **CWBcmapDBZ** </font> | Colormap, colors, contour levels and norm of CWB dbz pictures. |
| <font color="#a77864"> **CWBcmapRain** </font> | Colormap, colors, contour levels and norm of CWB accumulated daily/hourly rainfall pictures. |


******
pseudo_coord
------
[document](./documents/pseudo_coord.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **lonlat2xy** </font> | Convert to distance-based coordinate. |
| <font color="#a77864"> **xy2lonlat** </font> | Convert to longtitude-latitude coordinate. |


******
specialvar
------
[document](./documents/specialvar.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **uv2vrvt_rt** </font> | Calculate Vr (radial wind) and Vt (tangential wind) on radius-theta coordinate. |
| <font color="#a77864"> **uv2vrvt_xy** </font> | Calculate Vr (radial wind) and Vt (tangential wind) on cartesian coordinate. |
| <font color="#a77864"> **inertial_stability_xy** </font> | Calculate (cyclinic) inertial stability at x-y (longtitude-latitude) coordinate. |
| <font color="#a77864"> **inertial_stability_rt** </font> | Calculate (cyclinic) inertial stability at cylindrical (radius-theta) coordinate. |


******
temporary
------
[document](./documents/temporary.md) 

| Function / Class | Description |
| :--------------- | :---------- |
| <font color="#a77864"> **TemporaryObj** </font> | Using this instance to collecte some temporary variables. |


******