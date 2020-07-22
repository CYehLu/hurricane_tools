import os
import warnings
import requests
import zipfile
import pandas as pd


class DownloadWarning(Warning):
    def __str__(self):
        return 'Download new best track data'


class JMAbstParser:   
    """Parse JMA best track information."""
    
    def __init__(self, mode='lite'):       
        """
        Select return mode.
        
        mode: str, {'lite', 'full', 'txt'}. Default is `lite`.
            Mode `lite` would only remain `time`, `grade`, `center latitude`, `center longtitude`
            and `minimum sea level pressure` information.
            Mode `full` would keep all information in JMA best track document.
            Mode `txt` would return a string.
        """
        if mode not in ['lite', 'full', 'txt']:
            raise ValueError(f'Unavailable mode: {mode}')
        self._mode = mode
        
        # check if file exists
        if not os.path.isfile('./bst_all/bst_all.txt'):
            self._download()
        
    def _download(self):
        """download JMA best track data"""
        warnings.warn(DownloadWarning())
        
        url = 'https://www.jma.go.jp/jma/jma-eng/jma-center/rsmc-hp-pub-eg/Besttracks/bst_all.zip'
        resp = requests.get(url)
        
        if not os.path.isdir('./bst_all'):
            os.mkdir('./bst_all')
        
        with open('./bst_all/bst_all.zip', 'wb') as f:
            f.write(resp.content)

        with zipfile.ZipFile('./bst_all/bst_all.zip', 'r') as zipf:
            zipf.extractall('./bst_all/')
        
    def _check_year(self, allcontent, id_=None, year=None):
        """check the search year, and download new data if necessary"""
        yy2year = lambda yy: '19' + yy if '51' <= yy <= '91' else '20' + yy
        
        last_line = allcontent[-1]
        
        # if last_line is all blank or empty string
        if not last_line.strip():
            i = -1
            while not last_line.strip():
                last_line = allcontent[i]
                i -= -1
                
        last_yy = last_line[:2]
        last_year = yy2year(last_yy)
        
        if id_:
            yy = str(id_)[:2]   # e.g id_='1013', this is 2010 and yy='10'
            year = yy2year(yy)
            
        if int(year) > int(last_year):
            # download new data
            self._download()
            with open('./bst_all/bst_all.txt') as file:
                allcontent = file.readlines()
        
        return allcontent
    
    def _search_by_id(self, id_):
        """
        return the index of the header line, and the number of content lines
        """
        with open('./bst_all/bst_all.txt') as file:
            allcontent = file.readlines()
            
        allcontent = self._check_year(allcontent, id_=id_)
        
        for idx_line, line in enumerate(allcontent):
            if line.startswith('66666') and line.split()[1] == id_:
                header = line
                num_line = int(line.split()[2])
                return allcontent, idx_line, num_line
            
        # if find nothing
        raise ValueError(f'Can not find TC with id {id_}')
    
    def _search_by_name(self, name, year):
        """
        return the index of the header line, and the number of content lines
        """
        with open('./bst_all/bst_all.txt') as file:
            allcontent = file.readlines()
            
        allcontent = self._check_year(allcontent, year=year)
        
        for idx_line, line in enumerate(allcontent):
            if line.startswith('66666') and name.lower() in line.lower():
                nextline = allcontent[idx_line+1]
                yy = nextline[:2]
                if str(year)[2:] == yy:
                    header = line
                    num_line = int(line.split()[2])
                    return allcontent, idx_line, num_line
                else:
                    continue
        
        # if find nothing
        raise ValueError(f'Can not find `{name} ({year})`')
            
    def _parse_header(self, header):
        headcols = header.split()
        
        self.header = {
            'international_id': int(headcols[1]),
            'TC_id': int(headcols[3]),
            'flag': int(headcols[5]),
            'diff_time': int(headcols[6]),
            'name': headcols[7],
            'last_revision': int(headcols[8])
        }
        
        self.header_info = {
            'international_id': (
                '<International number ID>  \n'
                ' Last two digits of calendar year followed by 2-digit serial '
                'number ID of the storm of Tropical Storm (TS) intensity or '
                'greater'
            ),
            'TC_id': (
                '<Tropical cyclone number ID>  \n'
                ' Serial number ID of the storm of intensity with maximum '
                'sustained wind speed of 28 kt (near gale) or  greater'
            ),
            'flag': (
                '<Flag of the last data line> \n '
                ' 0 : Dissipation, \n'
                ' 1 : Going out of the responsible area of Japan Meteorological Agency (JMA)'
            ),
            'diff_time': (
                '<Difference between the time of the last data and the time of the final analysis>  \n'
                ' Unit : hour'
            ),
            'name': '<Name of the storm>',
            'last_revision': '<Date of the latest revision>'
        }
        
    def _parse_content(self, content):
        slice_dict = {
            'time': slice(0, 8),
            #'indicator': slice(9, 12),
            'grade': slice(13, 14), 
            'lat_center': slice(15, 18),
            'lon_center': slice(19, 24),
            'pres_center(hPa)': slice(24, 28),
            'max_sustained_ws(kt)': slice(33, 36),
            'dir_longest_radius_50kt_ws': slice(41, 42),
            'longest_radius_50kt_ws': slice(42, 46),
            'shortest_radius_50kt_ws': slice(47, 51),
            'dir_longest_radius_30kt_ws': slice(52, 53),
            'longest_radius_30kt_ws': slice(54, 57),
            'shortest_radius_30kt_ws': slice(58, 62),
            'landfall': slice(71, 72)
        }

        data = []
        for line in content:
            line_data = []
            for ss in slice_dict.values():
                line_data.append(line[ss])
            data.append(line_data)

        df = pd.DataFrame(data, columns=slice_dict.keys())
        
        self.content_info = {
            'time': '<Time of analysis> (UTC)',
            'grade': (
                '<Grade>  \n'
                ' 1 : Not used  \n'
                ' 2 : Tropical Depression (TD)  \n'
                ' 3 : Tropical Storm (TS)  \n'
                ' 4 : Severe Tropical Storm (STS)  \n'
                ' 5 : Typhoon (TY)  \n'
                ' 6 : Extra-tropical Cyclone (L)  \n'
                ' 7 : Just entering into the responsible area of JMA  \n'
                ' 8 : Not used  \n'
                ' 9 : Tropical Cyclone of TS intensity or higher  \n'
            ),
            'lat_center': '<Latitude of the center> Unit : degree',
            'lon_center': '<Longtitude of the center> Unit : degree',
            'pres_center(hPa)': '<Central pressure> Unit : hPa',
            'max_sustained_ws(kt)': '<Maximum sustained wind speed> Unit : knot (kt)',
            'dir_longest_radius_50kt_ws': (
                '<Direction of the longest radius of 50kt winds or greater>  \n'
                ' 1 : Northeast (NE)  \n'
                ' 2 : East (E)  \n'
                ' 3 : Southeast (SE)  \n'
                ' 4 : South (S)  \n'
                ' 5 : Southwest (SW)  \n'
                ' 6 : West (W)  \n'
                ' 7 : Northwest (NW)  \n'
                ' 8 : North (N)  \n'
                ' 9 : (symmetric circle)'
            ),
            'longest_radius_50kt_ws': '<The longest radius of 50kt winds or greater> Unit : nautical mile (nm)',
            'shortest_radius_50kt_ws': '<The shortest radius of 50kt winds or greater> Unit : nautical mile (nm)',
            'dir_longest_radius_30kt_ws': (
                '<Direction of the longest radius of 30kt winds or greater>  \n'
                ' 1 : Northeast (NE)  \n'
                ' 2 : East (E)  \n'
                ' 3 : Southeast (SE)  \n'
                ' 4 : South (S)  \n'
                ' 5 : Southwest (SW)  \n'
                ' 6 : West (W)  \n'
                ' 7 : Northwest (NW)  \n'
                ' 8 : North (N)  \n'
                ' 9 : (symmetric circle)'
            ),
            'longest_radius_30kt_ws': '<The longest radius of 30kt winds or greater> Unit : nautical mile (nm)',
            'shortest_radius_30kt_ws': '<The shortest radius of 30kt winds or greater> Unit : nautical mile (nm)',
            'landfall': (
                '<Indicator of landfall or passage>  \n'
                ' Landfall or passage over the Japanese islands occurred within '
                'one hour after the time of the analysis with this indicator.'
            )
        }
        
        return df
        
    def parse(self, id_=None, name=None, year=None):            
        """
        Parse JMA best track based on the "TC id" or "TC name + year".
        
        Parameters:
        ----------
        id_ : str
            International number ID
        name : str
            TC name.
        year : int
            The year the TC occurred.
        """
        if id_:
            allcontent, idx_line, num_line = self._search_by_id(id_)
        else:
            allcontent, idx_line, num_line = self._search_by_name(name, year)
        
        ### parse header line
        header = allcontent[idx_line]
        self._parse_header(header)
        
        ### parse content
        content = allcontent[idx_line+1: idx_line+1+num_line]
        
        if self._mode == 'txt':
            self.content = content
            return content
        
        else:
            df = self._parse_content(content)
            
            # check if TC is between 1951~1999 or 2000~, and convert to dt
            if df['time'].between('51', '99').all():
                df['time'] = '19' + df['time']
            else:
                df['time'] = '20' + df['time']
            df['time'] = pd.to_datetime(df['time'], format='%Y%m%d%H')
            
            # convert str to int and convert lat / lon
            df = df.applymap(lambda s: int(s) if isinstance(s, str) and s.isdigit() else s)
            df['lat_center'] = df['lat_center'].astype(float) * 0.1
            df['lon_center'] = df['lon_center'].astype(float) * 0.1

            if self._mode == 'lite':
                df = df.iloc[:,:5]
                
            self.content = df
            return df