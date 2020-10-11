import numpy as np


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