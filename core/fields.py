
from math import sqrt
from pyMD.misc import sampler

class field(object):
    
    def addRegion(self, region):
        self.region = region
        
    def update(self):
        raise NotImplementedError("All fields must have an update function to work.")
    
        
#==============================================================================
#         Thermostats
#==============================================================================
        

class thermostat(field):       

    
    def __init__(self, bathT, dt):
        self.bathT = bathT
        self.dt = dt
        self.T = None 

    def apply(self, atom):
        raise NotImplementedError("All thermostats must have an apply function.")

    def updateLocal(self):
        pass

    def update(self):
        
        self.T = sampler.getTemperature(self.region.atoms)
        self.updateLocal()  
    
        print len(self.region.atoms), sampler.getNumberOfFrees(self.region.atoms)    
    
        for atom in self.region.atoms:
            if atom.sticky:
                continue
            
            self.apply(atom)

class byScaleVelocity(thermostat):
    
    gamma = None

    def __init__(self, bathT, dt, tau):
        
        self.tau = tau
        
        super(byScaleVelocity, self).__init__(bathT, dt)
        

    def apply(self, atom):
         atom.vel *= self.gamma
    
    def updateLocal(self):
        self.gamma = self.getGamma()

    def getGamma(self):
        raise NotImplementedError("All velocity scaling thermostats must implement a method for calculating the scale.")


class BerendsenThemostat(byScaleVelocity):
    
    def getGamma(self):
        
        return sqrt(1 + self.dt/self.tau*(self.bathT/self.T - 1))
        
    
#==============================================================================
#       Pressure
#==============================================================================


class pressureField(field):

    def update(self):
        
        self.P = sampler.getPressure(self.region.atoms, V=self.region.volume)
        print self.region, self.P



