# temporary  

[[source](.././hurricane_tools//temporary.py)]  

class <span style="color:#a77864">**TemporaryObj**</span>

    Using this instance to collecte some temporary variables.


| Methods | Description |
| :------ | :---------- |
| <font color="#a77864"> **\_\_init\_\_** </font> | Set attributes and corresponding values. |
| <font color="#a77864"> **add** </font> | Add attribute. |
| <font color="#a77864"> **get** </font> | Get the value of attribute by its name. |
| <font color="#a77864"> **delete** </font> | Delete attribute by its name. |
| <font color="#a77864"> **close** </font> | delete all attribute. |


<span style="color:#cca99b">TemporaryObj</span>.<span style="color:#a77864">**\_\_init\_\_**</span>**(self, \*\*kwargs)**

        Set attributes and corresponding values.

  
<span style="color:#cca99b">TemporaryObj</span>.<span style="color:#a77864">**add**</span>**(self, name, val)**

        Add attribute. 
        
        These two lines are equivalent:
        >>> tmp = TemporaryObj(); tmp.add('a', 10) 
        >>> tmp = TemporaryObj(); tmp.a = 10

  
<span style="color:#cca99b">TemporaryObj</span>.<span style="color:#a77864">**get**</span>**(self, name)**

        Get the value of attribute by its name.
        
        These two lines are equivalent:
        >>> tmp = TemporaryObj(a=10); print(tmp.get('a'))
        >>> tmp = TemporaryObj(a=10); print(tmp.a)

  
<span style="color:#cca99b">TemporaryObj</span>.<span style="color:#a77864">**delete**</span>**(self, name)**

        Delete attribute by its name.

  
<span style="color:#cca99b">TemporaryObj</span>.<span style="color:#a77864">**close**</span>**(self)**

        delete all attribute.

  
******