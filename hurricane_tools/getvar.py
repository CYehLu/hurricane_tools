import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset, MFDataset

_HAS_IMPORT_F90GETVAR = False

try:
    from .fortran.f90tk import calc_tk
    from .fortran.f90slp import dcomputeseaprs, dcomputeseaprs_nt
    from .fortran.f90dbz import calcdbz, calcdbz_nt
    from .fortran.f90vort import dcomputeabsvort, dcomputepv, dcomputeabsvort_nt, dcomputepv_nt
    _HAS_IMPORT_F90GETVAR = True
except ModuleNotFoundError:
    pass


__all__ = [
    'destagger',
    'GetVar'
]


def destagger(var, stagger_dim):
    """
    Convert variable from staggered to unstagger grid.
    
    This is almost exactly the same as `wrf.destagger`, but I simplified to 
    destagger process and restricted the `stagger_dim` to one of z, y, or, x
    dimension to improve the efficiency.
    
    Parameter
    ---------
    var : ndarray
        The variable on the staggered grid. The dimensions of this variable
        must be (..., z, y, x) or (..., y, x). 
    stagger_dim : int
        The dimension index to destagger.
        Only the rightmost 3 dimensions (z, y, or x) can be destaggered.
    
    Return
    ------
    Variable on the unstaggered grid.
    """
    ndim = var.ndim
    
    # sdim: -1 (x stagger), -2 (y stagger) or -3 (z stagger)
    sdim = stagger_dim - ndim if stagger_dim >= 0 else stagger_dim
        
    if sdim == -1:
        return (var[...,1:] + var[...,:-1]) / 2
    elif sdim == -2:
        return (var[...,1:,:] + var[...,:-1,:]) / 2
    elif sdim == -3:
        return (var[...,1:,:,:] + var[...,:-1,:,:]) / 2
    else:
        raise ValueError("Unavailable stagger_dim (must be one of the three rightmost dimensions).")


class GetVar:
    """
    Get variables. It is similar to `wrf.getvar` by wrf-python, but here I 
    store every intermediate variables, reduce the amount of function calling
    and rewrite the fortran functions to speed up.
    """
    
    def __init__(self, filename, timeidx=None):
        """
        Initialization the instance.

        Parameters
        ----------
        filename : str
            The netCDF file name
        timeidx : optional, int, list, or slice()
            The time index of variables.
            Default is None, and it would use `slice(0, None)`, the all time index.
        """
        if not _HAS_IMPORT_F90GETVAR:
            raise ModuleNotFoundError("The Fortran modules should be compiled before using it.")

        self.variables = {}
        self.filename = filename
        self.ncfile = Dataset(filename)

        if timeidx is None:
            self.timeidx = slice(0, None)
        elif isinstance(timeidx, (slice, int, list)):
            self.timeidx = timeidx
        else:
            raise ValueError(f"Unavailable timeidx type: {type(timeidx)}")
        
        # define f90 func for diagnosis variables
        if self.ncfile.dimensions['Time'].size == 1:
            self._func = {
                'tk': calc_tk,
                'slp': dcomputeseaprs,
                'dbz': calcdbz,
                'avo': dcomputeabsvort,
                'pvo': dcomputepv
            }
            
        else:
            self._func = {
                'tk': calc_tk,     
                'slp': dcomputeseaprs_nt,
                'dbz': calcdbz_nt,
                'avo': dcomputeabsvort_nt,
                'pvo': dcomputepv_nt
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
        
        if variables:
            del self.variables
            
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        self.close(ncfile=True, variables=True)
        
    def _tk(self, func_tk):
        p = self.get('P')
        pb = self.get('PB')
        pres = np.squeeze(p + pb)
        
        t = self.get('T')
        theta = np.squeeze(t + 300)
        
        # ravel to 1-d array and convert to fortran type
        shape = pres.shape
        pres = np.asfortranarray(pres.ravel())
        theta = np.asfortranarray(theta.ravel())
        
        # f90tk.calc_tk only allow 1-d array (for the efficiency)
        tk = func_tk(pres, theta).reshape(shape)
        return tk

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
        ph = destagger(ph, -3)

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
    
    def _pres(self):
        return 0.01 * (self.get('P') + self.get('PB'))
    
    def _geopt(self):
        return self.get('PH') + self.get('PHB')
    
    def _theta(self):
        return self.get('T') + 300
    
    def _dbz(self, func_dbz, use_varint=False, use_liqskin=False):
        # get necessary variables
        pres = np.squeeze(self.get('pres')) * 100   # hPa to Pa
        tk = np.squeeze(self.get('tk'))
        qv = np.squeeze(self.get('QVAPOR'))
        qr = np.squeeze(self.get('QRAIN'))
        
        try:
            qs = np.squeeze(self.get('QSNOW'))
        except KeyError:
            qs = np.zeros(qv.shape, dtype=qv.dtype)
            
        try:
            qg = np.squeeze(self.get('QGRAUP'))
        except KeyError:
            qg = np.zeros(qv.shape, dtype=qv.dtype)
            
        # convert to fortran type and shape from 'zyx' to 'xyz'
        pres_f = np.asfortranarray(pres.T)
        tk_f = np.asfortranarray(tk.T)
        qv_f = np.asfortranarray(qv.T, dtype=np.float64)
        qr_f = np.asfortranarray(qr.T, dtype=np.float64)
        qs_f = np.asfortranarray(qs.T, dtype=np.float64)
        qg_f = np.asfortranarray(qg.T, dtype=np.float64)
        
        # `sn0` = 0 if qs is all 0
        sn0 = 1 if qs.any() else 0
        
        ivarint = 1 if use_varint else 0
        iliqskin = 1 if use_liqskin else 0
        
        dbz_f = func_dbz(pres_f, tk_f, qv_f, qr_f, qs_f, qg_f, sn0, ivarint, iliqskin)
        dbz = np.asanyarray(dbz_f.T, order='c')
        return dbz
    
    def _avo(self, func_avo):
        # read variables
        u = np.squeeze(self.get('U'))
        v = np.squeeze(self.get('V'))
        msfu = np.squeeze(self.get('MAPFAC_U'))
        msfv = np.squeeze(self.get('MAPFAC_V'))
        msfm = np.squeeze(self.get('MAPFAC_M'))
        f = np.squeeze(self.get('F'))
        
        dx = self.ncfile.DX
        dy = self.ncfile.DY
        
        # convert to fortran type and shape from 'zyx' to 'xyz'
        u_f = np.asfortranarray(u.T)
        v_f = np.asfortranarray(v.T)
        msfu_f = np.asfortranarray(msfu.T)
        msfv_f = np.asfortranarray(msfv.T)
        msfm_f = np.asfortranarray(msfm.T)
        f_f = np.asfortranarray(f.T)
        
        avo_f = func_avo(u_f, v_f, msfu_f, msfv_f, msfm_f, f_f, dx, dy)
        avo = np.asanyarray(avo_f.T, order='C')
        return avo
    
    def _pvo(self, func_pvo):
        u = np.squeeze(self.get('U'))
        v = np.squeeze(self.get('V'))
        theta = np.squeeze(self.get('T')) + 300
        pres = np.squeeze(self.get('pres')) * 100   # hPa to Pa
        msfu = np.squeeze(self.get('MAPFAC_U'))
        msfv = np.squeeze(self.get('MAPFAC_V'))
        msfm = np.squeeze(self.get('MAPFAC_M'))
        f = np.squeeze(self.get('F'))
        
        dx = self.ncfile.DX
        dy = self.ncfile.DY
        
        # convert to fortran type and shape from 'zyx' to 'xyz'
        u_f = np.asfortranarray(u.T)
        v_f = np.asfortranarray(v.T)
        theta_f = np.asfortranarray(theta.T)
        pres_f = np.asfortranarray(pres.T)
        msfu_f = np.asfortranarray(msfu.T)
        msfv_f = np.asfortranarray(msfv.T)
        msfm_f = np.asfortranarray(msfm.T)
        f_f = np.asfortranarray(f.T)
        
        pvo_f = func_pvo(u_f, v_f, theta_f, pres_f, msfu_f, msfv_f, msfm_f, f_f, dx, dy)
        pvo = np.asanyarray(pvo_f.T, order='C')
        return pvo
    
    def get_times(self, dt=None):
        """
        Get `Times` variable from netCDF file.
        
        Directly get `Times` from netCDF file would return a (ntime, 19) array, each column
        of array is a `numpy.bytes_` and it is not easy to read.
        This method would transform this np.bytes_ array, and return a list of string (or 
        datetime object, pandas DatetimeIndex) which is more convenient to read.
        
        Paramter
        --------
        dt : str, 'python' or 'pandas'. optional
            Default is None, it would return a list of string.
            If 'python', it would return a list of datetime object.
            If 'pandas', it would return a pandas.DatatimeIndex
        """
        times = self.get('Times')
        
        res = []
        for tt in times:
            res.append(''.join(list(map(lambda b: b.decode('utf-8'), tt))))

        if dt is None:
            return res

        elif dt == 'python':
            return list(map(lambda s: datetime.datetime.strptime(s, '%Y-%m-%d_%H:%M:%S'), res))

        elif dt == 'pandas':
            return pd.to_datetime(res, format='%Y-%m-%d_%H:%M:%S')

        else:
            raise ValueError(f"Argument `dt` should be None, 'python', or 'pandas'. Not `{dt}`.")
        
    def get(self, var_name, squeeze=True, filled=True):
        """
        Get variable by its name.

        Parameters
        ----------
        var_name : str
            Variable name.
            It can be the variable name of the netCDF file, or some diagnosis variables
            list below:
                'slp'   --  Sea Level Pressure
                'tk'    --  Temperature (unit: K)
                'pres'  --  Pressure (unit: hPa)
                'geopt' --  Geopotential (unit: m2 s-2)
                'theta' --  (Full) Potential Temperature (unit: K)
                'dbz'   --  Radar Reflectivity 
                'avo'   --  Absolute Vorticity (unit: 10-5 s-1)
                'pvo'   --  Potential Vorticity (unit: PVU)
                
        squeeze : bool, default is True
            Determine whether to squeeze the returned variable.
            
        filled : bool or scalar. Default is True.
            Determine the filled value of returned variable, if it is a `MaskArray`.
            
            If filled = True:
                Filled np.nan to variable. Returned variable type is `ndarray`.
            If filled = False:
                Do not fill any value to variable.
                Returned variable type is `MaskArray` (if wrf-output-variable) or `ndarray` (if 
                diagnois variable).
            If filled = scalar (int or float):
                Filled given value to variable. Returned variable type is `ndarray`.
                When filled = np.nan, it is equivalent to filled = True.
        """
        
        if var_name in self.variables.keys():
            # get variable from cache
            var = self.variables[var_name]
        
        else:
            if var_name in self.ncfile.variables.keys():
                # read variable from (wrfout) netCDF file
                var = self.ncfile.variables[var_name][self.timeidx]
                
            else:
                # calculate diagnosis variable
                if var_name == 'slp':
                    var = self._slp(self._func['slp'])
                elif var_name == 'tk':
                    var = self._tk(self._func['tk'])
                elif var_name == 'pres':
                    var = self._pres()
                elif var_name == 'geopt':
                    var = self._geopt()
                elif var_name == 'theta':
                    var = self._theta()
                elif var_name == 'dbz':
                    var = self._dbz(self._func['dbz'])
                elif var_name == 'avo':
                    var = self._avo(self._func['avo'])
                elif var_name == 'pvo':
                    var = self._pvo(self._func['pvo'])
                else:
                    raise KeyError(f"Unavailable variable: {var_name}")
            
            # update cache
            self.variables[var_name] = var
            
            
        # post-process returned variable
        if squeeze:
            var = np.squeeze(var)
            
        if isinstance(var, np.ma.MaskedArray):
            if isinstance(filled, bool):
                if filled:
                    var = var.filled(np.nan)
            elif isinstance(filled, (int, float)):
                var = var.filled(filled)
            else:
                raise ValueError(f"`filled` can only be bool, int or float, not {type(filled)}")
        
        return var