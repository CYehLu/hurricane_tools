#!/bin/sh

echo -n "Delete .so files... "
rm f90*.so
echo "[done]"

# for class GetVar
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp slp.f90 -m f90slp
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp tk.f90 -m f90tk
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp dbz.f90 -m f90dbz
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp vort.f90 -m f90vort

# for class Interpz3d
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp interp.f90 -m f90interp
