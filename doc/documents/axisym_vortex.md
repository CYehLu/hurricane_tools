# axisym_vortex  

[[source](../.././hurricane_tools//axisym_vortex.py)]  

<span style="color:#a77864">**rankine_vortex**</span>**(r, vmax, rmax, alpha=1)**

    Classic rankine vortex.
    
    V(r) = |- Vmax * (r / Rmax) ,       if r <= Rmax
           |- Vmax * (r / Rmax)**(-£\) , if r > Rmax
    in the classic case, £\ = 1.
    
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



******
<span style="color:#a77864">**holland80**</span>**(r, pc, pn, A, B, rho=None, f=None)**

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
        Mon. Wea. Rev., 108, 1212¡V1218
        https://reurl.cc/D9rbZd



******
<span style="color:#a77864">**willoughby06**</span>**(r, vmax, rmw, trans_width, n, decay_length1, decay_length2=None, proportion=None)**

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



******