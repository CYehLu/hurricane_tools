import numpy as np
from scipy.optimize import root


__all__ = [
    'rankine_vortex',
    'holland80',
    'willoughby06'
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


def _rhs(r1, rmw=None, transition_width=None):
    xi = (rmw - r1) / transition_width
    if xi <= 0:
        return 0
    elif xi >= 1:
        return 1
    else:
        return 126*xi**5 - 420*xi**6 + 540*xi**7 - 315*xi**8 + 70*xi**9
    
def _vt(r, vmax, rmw, r1, r2, n, x1, x2, a):
    v_inner = lambda r: vmax * (r/rmw)**n
    v_outer = lambda r: vmax * ((1-a)*np.exp(-(r-rmw)/x1) + a*np.exp(-(r-rmw)/x2))

    if isinstance(r, (int, float)):
        r = np.array([r])
    elif isinstance(r, (list, tuple)):
        r = np.array(r)
    elif isinstance(r, np.ndarray):
        pass

    inner_idx = (r >= 0) & (r <= r1)
    trans_idx = (r > r1) & (r < r2)
    outer_idx = r >= r2

    result = np.empty_like(r)
    result[inner_idx] = v_inner(r[inner_idx])
    result[outer_idx] = v_outer(r[outer_idx])

    xi = (r - r1) / (r2 - r1)
    w = 126*xi**5 - 420*xi**6 + 540*xi**7 - 315*xi**8 + 70*xi**9
    result[trans_idx] = (1-w[trans_idx]) * v_inner(r[trans_idx]) + w[trans_idx] * v_outer(r[trans_idx])
    
    return result 

def willoughby06(r, vmax, rmw, trans_width, n, decay_length1, decay_length2=None, proportion=None):
    """
    The empirical formula of TC axisymmetric structure by Willoughby (2006).

    Parameters
    ----------
    r :  scalar, 1-d array-like or 2-d array-like
        Radius to TC center
    vmax : scalar
        Maximum tangential wind speed
    rmw : scalar
        Radius of maximum wind
    trans_width : scalar
        Transition width between inner and outer vortex
    n : scalar
        Exponent for the power law inside the eye
    decay_length1 : scalar
        The first exponential decay length in the outer vortex
    decay_length2 : scalar, optional
        The second exponential decay length in the outer vortex
    proportion : scalar, optional
        The mix proportion of decay_length1 and decay_length in the outer vortex

    Return
    ------
    Vt : scalar, 1-d array-like or 2-d array-like
        Tangential wind speed.

    Note
    ----
    If `proportion` is 0, it would decay to single-exponential profile.
    
    Reference
    ---------
    [1] H. E. Willoughby, R. W. R. Darling, M. E. Rahn, "Parametric Representation of the Primary 
        Hurricane Vortex. Part II: A New Family of Sectionally Continuous Profiles"
        Mon. Wea. Rev., 134, 1102-1120
        https://journals.ametsoc.org/view/journals/mwre/134/4/mwr3106.1.xml
    """
    if decay_length2 is None:
        decay_length2 = 1
    if proportion is None:
        proportion = 0
    
    # symbols used in the paper
    a = proportion
    x1 = decay_length1
    x2 = decay_length2

    lhs = (n*(1-a)*x1 + n*a*x2) / (n*(1-a)*x1 + n*a*x2 + rmw)
    f = lambda r1: lhs - _rhs(r1, rmw, trans_width)
    res = root(f, x0=rmw-trans_width/2)  # the initial value is chosen to make xi(r=rmw) be 1/2

    if not res.success:
        raise ValueError("Solve Fail")
    
    r1 = res.x
    r2 = trans_width + r1
    v = _vt(r, vmax, rmw, r1, r2, n, x1, x2, a)
    return v 