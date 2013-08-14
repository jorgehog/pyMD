
from numpy import empty, random as npRandom, array, zeros

from random import uniform

from pyMD.misc.parser import tableParser

import sys


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
        self.EkPrev = None
        
        
    def initialize(self, mixingProperties, mesh):
        
        self.mesh = mesh
        mixTable = self.unpackMix(mixingProperties)  
      
        for i, _atom in enumerate(self.atoms):
            _atom.initialize(mesh, *mixTable[i])
            
        self.cancelLinearMomentum()
        
        self.getKineticEnergy()       
        
        self.calculateForces()
        
    def setSimulator(self, sim):
        self.simulator = sim
    
    def unpackMix(self, mixProp):
        
        unpacker = tableParser(self, self.mesh, self.N)
        
        table = unpacker.unpackMix(mixProp)
        
        return table
    
    def getColor(self, atom):
        
        if atom.sticky:
            sigmas = self.fixedSigmas
        else:
            sigmas = self.freeSigmas
            
        return sigmas.index(atom.sigma)        
        
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
#                print "r_rel", relPos
#                print "v", atom1.vel, atom2.vel
                
                force = self.forceModel.calculateForce(atom1, atom2, relPos, relPos2)
                
                atom1.updateForce(force)
                atom2.updateForce(-force)
         
        self.checkForces()
  
       
    def getTotalLinearMomentum(self):
        
        pTot = zeros(shape=self.mesh.dim)
        
        N = 0
        for _atom in self.atoms:
            if _atom.sticky:
                continue
            
            pTot += _atom.mass*_atom.vel
            N += 1
        
        return pTot, N
    
    
    def cancelLinearMomentum(self):
        
        pTotScaled, N = self.getTotalLinearMomentum()
        
        for _atom in self.atoms:
            if _atom.sticky:
                continue
            
                _atom.vel -= pTotScaled/(_atom.mass*N)
        
        self.checkLinearMomentum()
 
        
    def checkForces(self):
        S = empty(shape=2,)
        S.fill(0)

        for atom in self.atoms:
            S += atom.force
        
        for ei in S:
            if ei > 1E-8:
                print "Round-off errors in forcesum. %g Breaking simulation" % ei
  
                self.simulator.stopped = True

                
    def checkLinearMomentum(self):
        
        pTot = self.getTotalLinearMomentum()[0].sum() 
        
        if pTot > 1E-10:
            print "Total linear momentum is not zero. %g Breaking simulation" % pTot

                
    def getKineticEnergy(self):
        
        Ek = 0    
        for _atom in self.atoms:
            if _atom.sticky:
                continue
            
            Ek += 0.5*_atom.mass*(_atom.vel**2).sum()
            
        self.checkKineticEnergy(Ek)
                
        return Ek

    
    def checkKineticEnergy(self, Ek):
        
        if not self.EkPrev:
            self.EkPrev = Ek
                
        else:
            if abs(Ek - self.EkPrev) > 1e-3:
                print "Energy not conserved. %g / %g Breaking simulation" % (self.EkPrev, Ek)
                
                self.simulator.stopped = True
                
            else:
                self.EkPrev = Ek       
        