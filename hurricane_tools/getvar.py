import numpy as np
import wrf
from netCDF4 import Dataset, MFDataset

from .fortran.f90tk import calc_tk, calc_tk_nd
from .fortran.f90slp import dcomputeseaprs, dcomputeseaprs_nt
from .fortran.f90interp import find_level_1, find_level_n, interpz3d_1, interpz3d_n
from .fortran.f90dbz import calcdbz, calcdbz_nt
from .fortran.f90vort import dcomputeabsvort, dcomputepv, dcomputeabsvort_nt, dcomputepv_nt


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
        
        # define f90 func
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
                'tk': calc_tk_nd,
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
    
    def _pres(self):
        return 0.01 * (self.get('P') + self.get('PB'))
    
    def _geopt(self):
        return self.get('PH') + self.get('PHB')
    
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

    def get(self, var_name):
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
                'dbz'   --  Radar Reflectivity 
                'avo'   --  Absolute Vorticity (unit: 10-5 s-1)
                'pvo'   --  Potential Vorticity (unit: PVU)
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
                elif var_name == 'pres':
                    var = self._pres()
                elif var_name == 'geopt':
                    var = self._geopt()
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
            
            return var

        
class Interpz3d:
    """
    Interpolating variables on pressure coordinate.
    """
    
    def __init__(self, pres, level, missing_value=np.nan):
        """
        Initialize with pressure and levels.
        
        Parameter
        ---------
        pres : 3-d array, shape = (nz, ny, nx)
            pressure
        level : scalar or 1-d array with shape = (nlev,)
            interpolated pressure levels
        missing_value : scalar, optional
            Assign missing value. Default is `np.nan`
        """
        self.pres = pres
        self.level = level
        self.missing_value = missing_value
        
        if isinstance(level, (int, float)):
            find_level_func = find_level_1
            self._interpz3d_func = interpz3d_1
        elif isinstance(level, np.ndarray):
            find_level_func = find_level_n
            self._interpz3d_func = interpz3d_n
        else:
            raise ValueError(f"Unavailable `level` type : {type(level)}")
            
        # convert to fortran type and shape from `zyx` to `xyz`
        self._pres_f = np.asfortranarray(pres).T
        
        # find level index, shape = (nx, ny, nlev)
        self._lev_idx = find_level_func(self._pres_f, level)
    
    def interp(self, *var):
        """
        Interpolate
        
        Parameter
        ---------
        *var : 3-d array, shape = (nz, ny, nx)
            Interpolated variables. Their shapes should be the same.
            
        Return
        ------
        If len(var) == 1, it would return a interpolated variable with shape = (nlev, ny, nx)
        If len(var) > 1, it would return a list which elements are interpolated variable and
        shape = (nlev, ny, nx).
        """
        lev_idx = self._lev_idx
        pres_f = self._pres_f
        interpz3d_func = self._interpz3d_func
        level = self.level
        missing_value = self.missing_value
        
        if len(var) == 1:
            v_f = np.asfortranarray(var[0]).T
            v_interp = interpz3d_func(v_f, pres_f, level, lev_idx)   # (nx, ny)
            v_interp[lev_idx == 0] = np.nan
            return v_interp.T
        
        else:
            var_interp = []
        
            for v in var:
                v_f = np.asfortranarray(v).T
                v_interp = interpz3d_func(v_f, pres_f, level, lev_idx)   # (nx, ny, nlev)
                v_interp[lev_idx == 0] = missing_value
                var_interp.append(v_interp.T)
            
            return var_interp
            