# parse_wrf_center_mesg  

[[source](.././hurricane_tools//parse_wrf_center_mesg.py)]  

<span style="color:#a77864">**parse_wrf_rsl_error**</span>**(wrfrun_path)**

    Parse `rsl.error.0000` outputed from WRF to get the center information.
    If parse file successfully, it would create `center.txt` in the current
    folder.
    
    Parameters:
    ----------
    wrfrun_path: str
        The path of WRF run. e.g: wrfrun_path = '../WRFV3/run/'



******
<span style="color:#a77864">**find_centers_nearest_time**</span>**(filename, target_dates, formats=None)**

    Find the hurricane center data at specific times from the file
    which contains center information based on WRF rsl.error.0000.
    
    Parameters:
    ----------
    filename: 
        str, the file name of WRF hurricane centers file.
        The file is obtained by `get_wrf_ty_centers`.
    target_dates: 
        list of pandas datetime objects or list of str.
        The nearest datetime and its data of the center data in the 'filename'
        file would be find out.
    formats:
        str. The string format of target_dates if they are str. Option.
        
    Return:
    ------
    DataFrame



******