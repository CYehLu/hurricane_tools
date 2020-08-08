import numpy as np
import wrf
from netCDF4 import Dataset, MFDataset

from .fortran.f90tk import calc_tk, calc_tk_nd
from .fortran.f90slp import dcomputeseaprs, dcomputeseaprs_nt


class GetVar:
    """
    Get variables. It is similar to `wrf.getvar` by wrf-python,
    but here I store every intermediate variables to speed up.
    """
    
    def __init__(self, filename, timeidx=None):
        """
        Initialization the instance.

        Parameters
        ----------
        filename : str
            The netCDF file name
        timeidx : optional, int or slice()
            The time index of variables.
            Default is None, and it would use `slice(0, None)`, the all time index.
        """
        self.variables = {}
        self.filename = filename
        self.ncfile = Dataset(filename)

        if timeidx is None:
            self.timeidx = slice(0, None)
        elif isinstance(timeidx, (slice, int)):
            self.timeidx = timeidx
        else:
            raise ValueError(f"Unavailable timeidx type: {type(timeidx)}")
        
        if self.ncfile.dimensions['Time'].size == 1:
            self._func = {
                'tk': calc_tk,
                'slp': dcomputeseaprs
            }
            
        else:
            self._func = {
                'tk': calc_tk_nd,
                'slp': dcomputeseaprs_nt
            }

    def close(self, ncfile=True, variables=False):
        """
        Close attribute of GetVar instance.

        Parameters
        ----------
        ncfile : optional, bool
            Close GetVar.ncfile or not. Default is True.
        variables : optional, bool
            Delete GetVar.variables or not. Default is False.
        """
        if ncfile:
            self.ncfile.close()
        elif variables:
            del self.variables
        
    def _tk(self, func_tk):
        p = self.get('P')
        pb = self.get('PB')
        pres = np.squeeze(p + pb)
        
        t = self.get('T')
        theta = np.squeeze(t + 300)
        
        # convert to fortran type, and shape from (nz, ny, nx) to (nx, ny, nz)
        pres = np.asfortranarray(pres.T)
        theta = np.asfortranarray(theta.T)
        
        tk = func_tk(pres, theta)
        return np.asanyarray(tk.T, order='c')

    def _slp(self, func_slp):
        # read necessary variables
        p = self.get('P').copy()
        pb = self.get('PB')
        qvapor = self.get('QVAPOR').copy()
        ph = self.get('PH').copy()
        phb = self.get('PHB')
        tk = self.get('tk')
        
        # some preprocess of variables
        p += pb
        qvapor[qvapor < 0] = 0
        #ph = (ph + phb) / 9.81
        np.add(ph, phb, out=ph)
        np.divide(ph, 9.81, out=ph)
        ph = wrf.destagger(ph, -3)

        # convert to fortran type, shape from (nt, nz, ny, nx) to (nx, ny, nz)
        ph = np.asfortranarray(np.squeeze(ph).T)
        p = np.asfortranarray(np.squeeze(p).T)
        qvapor = np.asfortranarray(np.squeeze(qvapor).T)
        tk = np.asfortranarray(tk.T)

        # calculate sea level pressure
        #nx, ny = p.shape[:2]
        #slp = np.empty((nx, ny), np.float64, order='F')
        slp = func_slp(ph, tk, p, qvapor)
        slp = np.asanyarray(slp.T, order='c')
        return slp

    def get(self, var_name):
        """
        Get variable by its name.

        Parameters
        ----------
        var_name : str
            Variable name.
            It can be the variable name of the netCDF file, or some diagnosis variables
            list below:
                'slp'  --  Sea Level Pressure
                'tk'   --  Temperature (unit: K)
        """
        
        if var_name in self.variables.keys():
            # get variable from cache
            var = self.variables[var_name]
            return var
        
        else:
            if var_name in self.ncfile.variables.keys():
                # read variable from (wrfout) netCDF file
                var = self.ncfile.variables[var_name][self.timeidx]
                
            else:
                # calculate diagnois variable from fortran
                if var_name == 'slp':
                    var = self._slp(self._func['slp'])
                elif var_name == 'tk':
                    var = self._tk(self._func['tk'])
                else:
                    raise ValueError(f"Unavailable variable: {var_name}")
            
            # update cache
            self.variables[var_name] = var
            
            return var
