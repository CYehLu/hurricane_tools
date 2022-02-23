# hurricane_tools
-----

Some calculation routines used in my tropical cyclone research.   
For example:
* calculate axisymmetric wind field
* transform data from cartesian to polar coordinate
* interpolate data from model eta vertical level to pressure level
* find radaius of maximum wind
* calculate intertial stability
* fourier decomposition
* ...etc

## Environment and dependencies

I only test it in Linux environment.  
This package depends on the following packages:
* numpy (1.19.5)
* scipy (1.4.1)
* matplotlib (3.3.1)
* pandas (1.0.3)
* netCDF4 (1.5.3)
* requests (2.14.2)  

The versions in the parentheses are which I developed and tested. It may works well under different versions if these packages are not changed drastically.


## Installation

It can be downloaded by

    $ git clone https://github.com/CYehLu/hurricane_tools.git
    
and installed by pip

    $ pip install -e ./hurricane_tools
    

### Compile Fortran files
Some functions will use Fortran subroutines defined in `hurricane_tools/fortran/`. The functions that will use Fortran subroutines are listed below:

* `hurricane_tools.coord_transform.XY2RT(..., intp='fortran')`
* `hurricane_tools.coord_transform.RT2XY(..., intp='fortran')`
* `hurricane_tools.coord_transform.Interpz3d`
* `hurricane_tools.getvar.GetVar`
* Any function with argument `intp='fortran'` or the global variable `hurricane_tools.config.INTERP_OPTION = 'fortran'`  

It is necessary to compile Fortran files before using these functions.  
It can be compiled by the following steps:

    $ cd hurricane_tools/hurricane_tools/fortran
    $ make
    
It will generate seven `*.so` files if the compilation is successful. Note that I only test it with `gfortran` compiler.  


## Documents
See the list of module and the corresponding functions and classes [here](./doc/table.md).  

The document of each module:  
* [axisym_vortex](./doc/documents/axisym_vortex.md)  
* [bst_parser](./doc/documents/bst_parser.md)  
* [center](./doc/documents/center.md)  
* [circule](./doc/documents/circular.md)  
* [coord_transform](./doc/documents/coord_transform.md)  
* [distance](./doc/documents/distance.md)  
* [getvar](./doc/documents/getvar.md)  
* [interpolate](./doc/documents/interpolate.md)  
* [plot](./doc/documents/plot.md)  
* [pseudo_coord](./doc/documents/pseudo_coord.md)  
* [specialvar](./doc/documents/specialvar.md)  
* [temporary](./doc/documents/temporary.md)  