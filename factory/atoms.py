
from numpy import empty, zeros, dot

from pyMD.misc.parser import tableParser
from pyMD.misc import sampler

class atom:
        
    def __init__(self):
        
        self.initialized = False
        
    def initialize(self, mesh, N, sticky, sigma, eps, mass, initPos, initVel):
        
        self.sigma, self.eps, self.mass = sigma, eps, mass
        
        self.sticky = sticky
        
        self.mesh = mesh
    
        self.pos = initPos
        self.vel = initVel
        self.force = empty(mesh.dim)
        self.virial = zeros(N)
        
        self.initialized = True
         

    def incrementVel(self, dV):
        
        self.vel += dV
    
    def incrementPos(self, dR):
        
        self.pos = self.mesh.checkNewPos(self.pos + dR, self.vel)
        
    def updateForce(self, newForce):
        
        self.force += newForce
        
    def resetForce(self):
        
        self.force.fill(0)
        


   

class ensemble:
    
    def __init__(self, N, forceModel, T):
        
        self.atoms = [atom() for i in range(N)]
        self.forceModel = forceModel
        self.N = N
        self.EkPrev = None
        self.T = T
        self.pTot = None
        
        
    def initialize(self, mixingProperties, mesh):
        
        self.mesh = mesh
        
        unpacker = tableParser(self, mesh, self.N)
        
        unpacker.unpackMix(mixingProperties)
            
        self.cancelLinearMomentum()
        
        self.getKineticEnergy()       
        
        self.calculateForces()
        
    def setSimulator(self, sim):
        self.simulator = sim
    
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
        
        self.virial = 0
        for i, atom1 in enumerate(self.atoms[:-1]):    
    
            for j, atom2 in enumerate(self.atoms[(i+1):]):
                if atom1.sticky and atom2.sticky:
                    continue
           
                relPos, relPos2 = self.getRelPos(atom1, atom2)
#                print "r_rel", relPos
#                print "v", atom1.vel, atom2.vel
                
                force = self.forceModel.calculateForce(atom1, atom2, relPos, relPos2)
                
                g_ij = dot(force, relPos)
                self.virial += g_ij

                atom1.virial[j] = g_ij
                atom2.virial[i] = g_ij                
                
                atom1.updateForce(force)
                atom2.updateForce(-force)
         
        self.virial/=(3.0*self.mesh.volume)
        
        self.checkForces()

    def getPressure(self):
        
        return sampler.getPressure(self.atoms, N=self.nFree, V=self.mesh.volume, T=self.T, virial=self.virial)
       
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
                print "Round-off errors in forcesum. %g Pausing simulation" % ei
#                raw_input()
#                self.simulator.stopped = True

                
    def checkLinearMomentum(self):
        
        pTot = self.getTotalLinearMomentum()[0].sum() 
        
        if pTot > 1E-10:
            print "Total linear momentum is not zero. %g Pausing simulation" % pTot
#            raw_input()
                
    def getKineticEnergy(self):
        
        Ek = sampler.getKineticEnergy(self.atoms)
            
        self.checkKineticEnergy(Ek)
        
    def getTemperature(self):

            self.T = sampler.getTemperature(self.atoms, N=self.nFree, Ek=self.EkPrev)
            
            
    
    def checkKineticEnergy(self, Ek):
        
        if not self.EkPrev:
            self.EkPrev = Ek
                
        else:
            if abs(Ek/self.EkPrev - 1) > 0.1:
                print "Energy not conserved. %g / %g Pausingsimulation" % (self.EkPrev, Ek)
#                raw_input()
#                self.simulator.stopped = True
                
            
            self.EkPrev = Ek       
        