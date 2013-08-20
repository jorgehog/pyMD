
from math import sqrt
from numpy import array
from pyMD.misc import sampler
from pyMD.factory.geometry import region as _reg

class field(object):
    
    _type = "field"
    _unit = ""
    
    region = _reg(0, 1, 0, 1)
    
    boundaryIndexMap = {"l" : 0, "r": 0, "t": 1, "b": 1}      
    normalVectors = {"l" : array([1,0]),
                              "r" : array([-1, 0]),
                              "t" : array([0, -1]),
                              "b" : array([0, 1])}
    
    def addRegion(self, region):
        self.region = region
        
    def update(self):
        raise NotImplementedError("All fields must have an update function to work.")
        
    def getMeasurement(self):
        return 0
        
    def __str__(self):
        
        s = [13, 14, 4]        
        
        return "<%s, %s, value: %.3f %s>" % (self._type.ljust(s[0]), 
                                           str(self.region).ljust(s[1]), 
                                        self.getMeasurement(), 
                                        self._unit.ljust(s[2]))
    
        
#==============================================================================
#         Thermostats
#==============================================================================
        

class thermostat(field):       

    _type = "thermostat"
    _unit = "K"
    
    def __init__(self, bathT, dt):
        self.bathT = bathT
        self.dt = dt
        self.T = 0 

    def apply(self, atom):
        raise NotImplementedError("All thermostats must have an apply function.")

    def updateLocal(self):
        pass

    def update(self):
        
        self.T = sampler.getTemperature(self.region.atoms)
        self.updateLocal()  
    
#        print len(self.region.atoms), sampler.getNumberOfFrees(self.region.atoms)    
    
        for atom in self.region.atoms:

            if atom.sticky:
                continue
            

            self.apply(atom)

    def getMeasurement(self):
        return self.T
        
        

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
        
        if self.T == 0:
            return 1
        
        return sqrt(1 + self.dt/self.tau*(self.bathT/self.T - 1))
        
    
#==============================================================================
#       Pressure
#==============================================================================


class pressureField(field):
    
    _type = "pressureField"
    _unit = "Pa"
    
    counter = 0

    def __init__(self, pExt=0, where=[], when=[-1, -1], dx=0.1, dy=0.1):
        
        self.dr = [dx, dy]
        
        self.pExt = pExt        
        
        if type(where) != list:
           where = [where]         
            
        self.where = where
        self.when = when
        
        
    def update(self):    
        
        if (self.when[0] <= self.counter <= self.when[1]):
#            for atom in self.region.atoms:
#                print atom.pos
#                
#            raw_input()
            for boundary in self.where:
                
                i = self.boundaryIndexMap[boundary]
                j = (i + 1)%2
                
                r = self.region.getShape()
                
                boundaryAtoms = self.region.extractBoundary(boundary, width=r[i]*self.dr[i])
                for atom in boundaryAtoms:
#                    print "atom at", atom.pos*1E9, "found in region", self.region.description, " at boundary", boundary, "with width", r[i]*self.dr[i]
#                    raw_input()
                    self.apply(atom, r[j], self.normalVectors[boundary])
#            print "----------------------------------"
            

        
        self.P = sampler.getPressure(self.region.atoms, V=self.region.volume)
        self.counter += 1

    def apply(self, atom, crossSectionArea, normalVector):
        
        f = self.pExt*crossSectionArea*normalVector
        
        print "adding", f, "to atom at", atom.pos
        atom.force += f
        atom.queuedForce += f
        


    def getMeasurement(self):
        return self.P
        



