

from math import sqrt
from numpy import array, zeros
from numpy.random import uniform, normal


class RNGs:
    
    finished = False


class randomOnMesh(RNGs):
    
        def __init__(self, mesh):
            
             self.shape = mesh.shape

        def __call__(self, **kwargs):
            
             plc = []
             
             for l in self.shape:
                 plc.append(uniform(0, l))
                 
             return array(plc)
             

class normalFromTermperature(RNGs):
    
    def __init__(self, T, dim = 2):

        self.dim = dim        
        self.T = T
        
    def __call__(self, **kwargs):
        
        return normal(scale = sqrt(self.T/kwargs["mass"]), size = self.dim)
        
        
class allZero(RNGs):
    
    def __init__(self, dim = 2):
        
        self.dim = dim
        
    def __call__(self, **kwargs):
        
        return zeros(self.dim)
        


