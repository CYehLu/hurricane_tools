__all__ = [
    'TemporaryObj'
]


class TemporaryObj:
    """Using this instance to collecte some temporary variables."""
    
    def __init__(self, **kwargs):
        """
        Set attributes and corresponding values.
        """
        for name, val in kwargs.items():
            self.__setattr__(name, val)
            
    def __repr__(self):
        return f'Temporary Object. Attributes : {list(self.__dict__.keys())}'
    
    def add(self, name, val):
        """
        Add attribute. 
        
        These two lines are equivalent:
        >>> tmp = TemporaryObj(); tmp.add('a', 10) 
        >>> tmp = TemporaryObj(); tmp.a = 10
        """
        self.__setattr__(name, val)
        
    def get(self, name):
        """
        Get the value of attribute by its name.
        
        These two lines are equivalent:
        >>> tmp = TemporaryObj(a=10); print(tmp.get('a'))
        >>> tmp = TemporaryObj(a=10); print(tmp.a)
        """
        return self.__getattribute__(name)
        
    def delete(self, name):
        """
        Delete attribute by its name.
        """
        self.__delattr__(name)
    
    def close(self):
        """
        delete all attribute.
        """
        # attrs: [('attr1', v1), ('attr2', v2), ...]
        attrs = list(self.__dict__.items())
        for name, val in attrs:
            self.__delattr__(name)