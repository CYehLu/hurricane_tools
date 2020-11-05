import os
import pandas as pd


__all__ = [
    'parse_wrf_rsl_error',
    'find_centers_nearest_time'
]


def parse_wrf_rsl_error(wrfrun_path):
    """
    Parse `rsl.error.0000` outputed from WRF to get the center information.
    If parse file successfully, it would create `center.txt` in the current
    folder.
    
    Parameters:
    ----------
    wrfrun_path: str
        The path of WRF run. e.g: wrfrun_path = '../WRFV3/run/'
    """
    command = 'grep ATCF ' + wrfrun_path + 'rsl.error.0000 > centers.txt'
    flag = os.system(command)

    if flag == 0:
        print('SUCCESSFUL. "centers.txt" is stored at the current folder.')
    else:
        print('FAIL')


def find_centers_nearest_time(filename, target_dates, formats=None):
    """
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
    """
    # read center file
    with open(filename) as cf:
        contents = cf.readlines()
        
    # only keep the time, lat and lon information
    contents = [line.split()[1:4] for line in contents]
    
    # parse contents
    df = pd.DataFrame(contents)
    df = df.set_index(0)
    df.index.name = None
    df.columns = ['lat', 'lon']
    df.index = pd.to_datetime(df.index, format='%Y-%m-%d_%H:%M:%S')
    df = df.applymap(lambda s: float(s))
    
    # convert to list of pd.datetime obj if "target_dates" are strings
    if isinstance(target_dates[0], str):
        target_dates = [pd.to_datetime(dt) for dt in target_dates]
        
    # find nearest dates
    idx = []
    for dt in target_dates:
        idx.append(df.index.get_loc(dt, method='nearest'))
        
    df2 = pd.DataFrame(
        data=df.iloc[idx,:].values,
        index=target_dates,
        columns=['lat', 'lon']
    )

    return df2

