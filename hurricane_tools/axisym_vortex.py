import numpy as np
from scipy.optimize import root


__all__ = [
    'rankine_vortex',
    'holland80',
    'willoughby04'
]


def rankine_vortex(r, vmax, rmax, alpha=1):
    """
    Classic rankine vortex.
    
    V(r) = |- Vmax * (r / Rmax) ,       if r <= Rmax
           |- Vmax * (r / Rmax)**(-α) , if r > Rmax
    in the classic case, α = 1.
    
    Parameters
    ----------
    r : scalar, 1-d array-like or 2-d array-like
        Radius
    vmax : scalar
        The maximum tangential wind speed.
    rmax : scalar
        The radius of maximum wind speed.
    alpha : scalar, optional
        The wind speed decrease rate outside Rmax. Default is 1.
        
    Return
    ------
    V : scalar, 1-d array-like or 2-d array-like
        Tangential wind speed.
    """
    if isinstance(r, (int, float)):
        r = np.array([r])
    elif isinstance(r, (list, tuple)):
        r = np.array(r)

    idx = r <= rmax
    v = np.empty_like(r)
    v[idx] = vmax * (r[idx] / rmax)
    v[~idx] = vmax * (r[~idx] / rmax)**(-alpha)
    return v


def holland80(r, pc, pn, A, B, rho=None, f=None):
    """
    The empirical formula of TC axisymmetric structure by Holland (1980).
    
    Follow this profile:
    radius of maximum wind = A ** (1/B)
    maximum wind speed = sqrt(B * (pn - pc) / (rho * e)), where e is nature logarithms.
    
    Parameters
    ----------
    r : scalar, 1-d array-like or 2-d array-like
        Radius
    pc : scalar
        The TC center pressure (Pa).
    pn : scalar
        The environmental pressure (Pa).
    A, B: scalar
        Scaling parameters (positive).
    rho : scalar, optional
        The air density. Default is 1.15 (kg/m^3).
    f : scalar, optional
        Coriolis coefficient. Default is 10**-5.
        
    Return
    ------
    V : scalar, 1-d array-like or 2-d array-like
        Tangential wind speed.
        
    Reference
    ---------
    [1] G. J. Holland, "An Analytic Model of the Wind and Pressure Profiles in Hurricanes"
        Mon. Wea. Rev., 108, 1212–1218
        https://reurl.cc/D9rbZd
    """
    if rho is None:
        rho = 1.15
    if f is None:
        f = 10 ** -5
    
    if isinstance(r, (int, float)):
        r = np.array([r])
    elif isinstance(r, (list, tuple)):
        r = np.array(r)
        
    v = np.zeros_like(r)
    not0 = ~np.isclose(r, 0)
    
    term1 = A * B * (pn - pc) * np.exp(-A / r[not0]**B) / (rho * r[not0]**B)
    term2 = r[not0] * f / 2
    v[not0] = np.sqrt(term1 + term2**2) - term2
    return v


def willoughby04(r, vmax, rmax, n, X1, X2, A, twidth):
    """
    The empirical formula of TC axisymmetric structure by Willoughby (2006).
    
    Reference
    ---------
    [1] H. E. Willoughby, R. W. R. Darling, M. E. Rahn, "Parametric Representation of the Primary 
        Hurricane Vortex. Part II: A New Family of Sectionally Continuous Profiles"
        Mon. Wea. Rev., 134, 1102-1120
        https://reurl.cc/Kk1EMM
    [2] Wu, C.-C., G.-Y. Lien, J.-H. Chen, and F. Zhang, "Assimilation of tropical cyclone track 
        and structure based on the Ensemble Kalman Filter (EnKF)"
        J. Atmos. Sci., 67, 3806-3822
        http://journals.ametsoc.org/doi/pdf/10.1175/2010JAS3444.1
    """
    raise ValueError('Unfinish')
    
    # determine `R1` and `R2`
    w_func = lambda xi: 126 * xi**5 - 480 * xi**6 + 540 * xi**7 - 315 * xi**8 + 70 * xi**9
    target = n * ((1-A) * X1 + A * X2) / (n * ((1-A) * X1 + A * X2) + rmax)
    target_func = lambda xi: w_func(xi) - target
    xi0 = root(target_func, 0.5).x
    R1 = rmax - xi0 * twidth
    R2 = R1 + twidth
    
    idx_i = r <= R1    # inner area
    idx_o = r >= R2    # outer area
    idx_t = (r > R1) & (r < R2)    # transition area
    
    vi_func = lambda r: vmax * (r / rmax)**n
    vo_func = lambda r: vmax * ( (1-A) * np.exp(-(r-rmax)/X1) + A * np.exp(-(r-rmax)/X2) )
    
    v = np.empty_like(r)
    v[idx_i] = vi_func(r[idx_i])
    v[idx_o] = vo_func(r[idx_o])
    
    xi = (r[idx_t] - R1) / (R2 - R1)
    w = w_func(xi)
    v[idx_t] = vi_func(r[idx_t]) * (1-w) + vo_func(r[idx_t]) * w
    
    return v