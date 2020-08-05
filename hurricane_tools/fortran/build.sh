#!/bin/sh

echo -n "Delete .so files... "
rm f90*.so
echo "[done]"
    
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp slp.f90 -m f90slp
f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp tk.f90 -m f90tk
