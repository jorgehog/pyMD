
from pyMD.core import simulator, forceModels, integrators

from pyMD.factory import atoms, geometry


shape = [2, 4]
periodicity = [False, True]
N = 100

mesh = geometry.rectMesh(periodicity, shape)

forceModel = forceModels.LennardJones()

atoms = atoms.ensemble(N, forceModel)

particleMixTable = {

                    "free" : {
                                "fraction"  : "remaining",
                                "nSpecies"  : 2,
                                "sigmas"    : [1, "rel 0 1e-5"],
                                "epses"     : [1, "as 0"],
                                "masses"    : [1, "as 0"],
                                
                                "positions" : "random uniform",
                                "velocities": "random normal"
                                
                    },
                             
                             
                    "fixed" : {
    
                                "fraction"  : "automatic",
                                
                                "nSpecies"  : "as free",                           
                                "sigmas"    : "as free",
                                "epses"     : "as free",  
                                "masses"    : "as free",
                                
                                "partition" : "cyclic",
                                
                                "positions" : [
                                               "boundary l", 
                                               "boundary r",
                                               "boundary t",
                                               "boundary b",
                                               "centerline x"
                                              ],
                                              
                                "thickness"  : 4,
                                "separation" : "rel ly 0.05"
                    
                    }

}
            

dt = 0.01
T = 10

integrator = integrators.EulerCramer(dt)

app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
app.run(T)







#TODO PARSER

#Start with fixed: Do positions, start with filling row by row and always check if occupied (brute force loop all fixed (within thickness of another) start with walls and loop iterate thickness in all other direction)



