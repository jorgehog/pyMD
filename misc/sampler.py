

                
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
        
