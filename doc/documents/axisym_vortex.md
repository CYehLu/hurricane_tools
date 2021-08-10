# axisym_vortex  

[[source](../.././hurricane_tools//axisym_vortex.py)]  

<span style="color:#a77864">**rankine_vortex**</span>**(r, vmax, rmax, alpha=1)**

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
        Mon. Wea. Rev., 108, 1212–1218
        https://reurl.cc/D9rbZd



******
<span style="color:#a77864">**willoughby04**</span>**(r, vmax, rmax, n, X1, X2, A, twidth)**

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



******