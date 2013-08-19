
from scipy.constants import k as kB
                
def getKineticEnergy(atoms):
    
    Ek = 0    
    for atom in atoms:
        if atom.sticky:
            continue
        
        Ek += 0.5*atom.mass*(atom.vel**2).sum()
        
    return Ek
    

def getNumberOfFrees(atoms):
    
    return len([atom for atom in atoms if not atom.sticky])
    

    
def getTemperature(atoms, Ek = None, N = None):
    
        if not Ek:
            Ek = getKineticEnergy(atoms)
            
        if not N:
            N = getNumberOfFrees(atoms)
    
        return 2*Ek/(3*N*1.38E-23)
        

def getVirial(atoms, V):
    
    N = len(atoms)    
    
    virial = 0
        
    for i in range(N):

        if atoms[i].sticky:
            continue
        
        for j in range(i+1, N):

            if atoms[j].sticky:
                continue

            virial += atoms[i].virials[j]
            
    virial /= 3.0*V
    
    return virial
    

def getPressure(atoms, V = 1, N = None, T=None, virial = None):
    
    if not N:
        N = getNumberOfFrees(atoms)

    #NB: Zero == False
    if T is None:
        T = getTemperature(atoms)
        
    if not virial:
        virial = getVirial(atoms, V)
        
    
    rho = N/float(V)
    
    return rho*kB*T + virial
    
    
    
    
    
    