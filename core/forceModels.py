
from math import sqrt


class forceModel():
    
    def calculateForce(self, atom1, atom2, relPos, relPos2):
        raise NotImplementedError("calculateForces method not implemented for class %s." \
                                % self.__class__.__name__)
                                

class LennardJones(forceModel):
    
    def calculateForce(self, atom1, atom2, relPos, relPos2):
        
        sigma = sqrt(atom1.sigma*atom2.sigma)
        eps = sqrt(atom2.eps*atom1.eps)
        
        l = sqrt(relPos2)
#        print l/sigma
        sigmaOverR6 = (sigma*sigma/relPos2)**3      
        
        return 4*eps*(sigmaOverR6**2 - sigmaOverR6)*relPos/l       
        

class LennardJonesReduced(forceModel):
    
    def calculateForce(self, atom1, atom2, relPos, relPos2):
        
        sigmaOverR6 = (1/relPos2)**3
        
        return 4*(sigmaOverR6**2 - sigmaOverR6)*relPos/sqrt(relPos2)
						

class CoulombLike(forceModel):
    
    def calculateForce(self, atom1, atom2, relPos, relPos2):
        
        eps = sqrt(atom2.eps*atom1.eps)
        
        return eps/relPos2
