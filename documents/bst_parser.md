# bst_parser  

[[source](.././hurricane_tools//bst_parser.py)]  

class <span style="color:#a77864">**DownloadWarning**</span>




| Methods | Description |
| :------ | :---------- |


******
class <span style="color:#a77864">**JMAbstParser**</span>

    Parse JMA best track information.


| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Select return mode. |
| <font color="#a77864"> **parse** </font> | Parse JMA best track based on the "TC id" or "TC name + year". |


<span style="color:#cca99b">JMAbstParser</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, mode='lite')**

        Select return mode.
        
        mode: str, {'lite', 'full', 'txt'}. Default is `lite`.
            Mode `lite` would only remain `time`, `grade`, `center latitude`, `center longtitude`
            and `minimum sea level pressure` information.
            Mode `full` would keep all information in JMA best track document.
            Mode `txt` would return a string.

  
<span style="color:#cca99b">JMAbstParser</span>.<span style="color:#a77864">**parse**</span>**(self, id_=None, name=None, year=None)**

        Parse JMA best track based on the "TC id" or "TC name + year".
        
        Parameters:
        ----------
        id_ : str
            International number ID
        name : str
            TC name.
        year : int
            The year the TC occurred.

  
******