

class integrator():
    
    def __init__(self, dt):
        self.dt = dt
        
    def getNewVel(self, atom):
        raise NotImplementedError("getNewVelShift is not implemented in class %s" % \
                                    self.__class__.__name__)
                                    
    def getNewPos(self, atom):
        raise NotImplementedError("getNewPosShift is not implemented in class %s" % \
                                    self.__class__.__name__)
                 
    def updateAtoms(self, ensemble):
        
        for atom in ensemble.atoms:

            if atom.sticky:
                continue 
    
            atom.incrementVel(self.getNewVelShift(atom))
            atom.incrementPos(self.getNewPosShift(atom))
            
        
        ensemble.calculateForces()
                                  

class EulerCramer(integrator):
    
    def getNewVelShift(self, atom):
        return atom.force / atom.mass * self.dt
        
    def getNewPosShift(self, atom):
        return atom.vel * self.dt
        
        
class VelocityVerlet(integrator):
    
    
    def updateAtoms(self, ensemble):
        
        for atom in ensemble.atoms:
            
            if atom.sticky:
                continue
            
            atom.incrementVel(self.getNewVelShift(atom))
            atom.incrementPos(self.getNewPosShift(atom))

        ensemble.calculateForces()  

        for atom in ensemble.atoms:
            
            if atom.sticky:
                continue
            
            atom.incrementVel(self.getNewVelShift(atom))
            
            
    
    def getNewVelShift(self, atom):
        return atom.force/(2*atom.mass)*self.dt 
    
    
    def getNewPosShift(self, atom):
        return atom.vel * self.dt