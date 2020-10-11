# destag  

[[source](.././hurricane_tools//destag.py)]  

<span style="color:#a77864">**destagger**</span>**(var, stagger_dim)**

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



******