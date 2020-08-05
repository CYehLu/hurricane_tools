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
    
    def __init__(self, filename):
        """Initialization. `filename` should be a str."""
        self.variables = {}
        self.filename = filename
        self.ncfile = Dataset(filename)
        
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

    def close(self):
        """close netCDF file instance"""
        self.ncfile.close()
        
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
        """Get variable by its name."""
        
        if var_name in self.variables.keys():
            # get variable from cache
            var = self.variables[var_name]
            return var
        
        else:
            if var_name in self.ncfile.variables.keys():
                # read variable from (wrfout) netCDF file
                var = self.ncfile.variables[var_name][:]
                
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