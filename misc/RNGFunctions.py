

from math import sqrt
from numpy import array, zeros
from numpy.random import uniform, normal


class RNGs:
    
    finished = False


class randomOnMesh(RNGs):
    
        def isOccupied(self, ens, pos, thresh):
            
            a = RNGs()
            setattr(a, "pos", pos)
            
            for atom in ens.atoms:
                if atom.initialized:
                    rel, rel2 = ens.getRelPos(atom, a)
                    
                    if rel2 < thresh*atom.sigma:
                        return True
            print "placed one!"
            return False

        def getNewPos(self, shape):
            
             plc = []
             
             for l in shape:
                 plc.append(uniform(0, l))
             
             return array(plc)          


        def __call__(self, **kwargs):
            
             ens = kwargs["ensemble"]
             
             pos = self.getNewPos(ens.mesh.shape)

             scale = 1

             k = 0
             while self.isOccupied(ens, pos, kwargs["sigma"]*scale):             
                 pos = self.getNewPos(ens.mesh.shape)
                 
                 if k > 10:
                     scale *= 0.9
                     k = 0
                 else:
                     k += 1
             
             return pos
             
             

class normalFromTermperature(RNGs):
    
    def __init__(self, T):

        self.T = T
        self.sqrtK = 3.7157E-12
            
    def __call__(self, **kwargs):
        
        return normal(scale = self.sqrtK*sqrt(self.T/kwargs["mass"]), size = kwargs["dim"])
        
        
class allZero(RNGs):
    
    def __init__(self, dim = 2):
        
        self.dim = dim
        
    def __call__(self, **kwargs):
        
        return zeros(self.dim)
        


