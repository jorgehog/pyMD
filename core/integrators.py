

class integrator():
    
    def __init__(self, dt):
        self.dt = dt
        
    def getNewVel(self, atom):
        raise NotImplementedError("getNewVel is not implemented in class %s" % \
                                    self.__class__.__name__)
                                    
    def getNewPos(self, atom):
        raise NotImplementedError("getNewVel is not implemented in class %s" % \
                                    self.__class__.__name__)
                 
    def updateAtoms(self, atoms):
        
        for atom in atoms:
            atom.incrementVel(self.getNewVelShift(atom))
            atom.incrementPos(self.getNewPosShift(atom))
                                    

class EulerCramer(integrator):
    
    def getNewVelShift(self, atom):
        return atom.force / atom.mass * self.dt
        
    def getNewPosShift(self, atom):
        return atom.vel * self.dt
        