
from matplotlib.pylab import scatter, show, clf
import time

class MDApp():
    
    def __init__(self, dt, mesh, atoms, mixTable, integrator):
        
        self.dt = dt
        self.mesh = mesh
        self.ensemble = atoms
        
        self.integrator = integrator       
        
        self.ensemble.initialize(mixTable, self.mesh)
        
    def run(self, T):
        

        n = int(T/self.dt)

      
 
        
        for i in range(n):
            
            self.ensemble.calculateForces()
            self.integrator.updateAtoms(self.ensemble.atoms)
            
            for atom in self.ensemble.atoms:
                scatter(atom.pos[0], atom.pos[1])
            show()
        
            time.sleep(0.01)
            clf()
            
                
            
            
        
        