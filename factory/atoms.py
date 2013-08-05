
from numpy import empty, random as npRandom, array

from random import uniform

from pyMD.misc.parser import tableParser


class atom():
        
    def initialize(self, mesh, sticky, sigma, eps, mass, initPos, initVel):
        
        self.sigma, self.eps, self.mass = sigma, eps, mass
        
        self.sticky = sticky
    
        self.pos = initPos
        self.vel = initVel
        self.force = empty(mesh.dim)
        
        self.mesh = mesh       

    def incrementVel(self, dV):
        
        self.vel += dV
    
    def incrementPos(self, dR):
        
        self.pos = self.mesh.checkNewPos(self.pos + dR)
        
    def updateForce(self, newForce):
        
        self.force += newForce
        
    def resetForce(self):
        
        self.force.fill(0)

   

class ensemble():
    
    def __init__(self, N, forceModel):
        
        self.atoms = [atom() for i in range(N)]
        self.forceModel = forceModel
        self.N = N
        
        
    def initialize(self, mixingProperties, mesh):
        
        self.mesh = mesh
        mixTable = self.unpackMix(mixingProperties)  
        
        
        for i, _atom in enumerate(self.atoms):
            _atom.initialize(mesh, *mixTable[i])
            
        
    
    def unpackMix(self, mixProp):
        
        unpacker = tableParser(self.mesh, self.N)
        
        table = unpacker.unpackMix(mixProp)
        
        return table
        
        
    def getRelPos(self, atom1, atom2):
        #SMARTIFY THIS
    
        minRelPos = atom1.pos - atom2.pos
        minRelPos2 = (minRelPos**2).sum()
        

     
        for i in range(self.mesh.dim):
            if self.mesh.p[i]:

                tmp = atom2.pos[i]
                atom2.pos[i] -= self.mesh.shape[i]
                
                test = ((atom1.pos - atom2.pos)**2).sum()
                if test < minRelPos2:
                    minRelPos = atom1.pos - atom2.pos
                    minRelPos2 = test
                    
                atom2.pos[i] += 2*self.mesh.shape[i]
                
                test = ((atom1.pos - atom2.pos)**2).sum()
                if test < minRelPos2:
                    minRelPos = atom1.pos - atom2.pos
                    minRelPos2 = test
                    
                atom2.pos[i] = tmp
                
        return minRelPos, minRelPos2
        
        
    def calculateForces(self):
        
        for atom in self.atoms:
            atom.resetForce()
        
        for i, atom1 in enumerate(self.atoms[:-1]):    
    
            for atom2 in self.atoms[(i+1):]:
                if atom1.sticky and atom2.sticky:
                    continue
           
                relPos, relPos2 = self.getRelPos(atom1, atom2)                
                
                force = self.forceModel.calculateForce(atom1, atom2, relPos, relPos2)
                
                atom1.updateForce(force)
                atom2.updateForce(-force)
        