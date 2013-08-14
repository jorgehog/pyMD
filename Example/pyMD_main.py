
from pyMD.core import simulator, forceModels, integrators

from pyMD.factory import atoms, geometry

from pyMD.misc.RNGFunctions import randomOnMesh, normalFromTermperature


shape = [1, 2]
periodicity = [False, True]
T = 1
N = 200

mesh = geometry.rectMesh(periodicity, shape)

forceModel = forceModels.LennardJones()
#forceModel = forceModels.CoulombLike()

atoms = atoms.ensemble(N, forceModel)

particleMixTable = {

                    "free" : {
                                "fraction"  : "remaining",
                                "nSpecies"  : 2,
                                "sigmas"    : [0.0001, "rel 0 0.99"],
                                "epses"     : [0.0001, "rel 0 0.99"],
                                "masses"    : [1, "rel 0 0.99"],
                                
                                "positions" : randomOnMesh(mesh),
                                "velocities": normalFromTermperature(T, len(shape))
                                
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
#                                               "boundary t",
#                                               "boundary b",
#                                               "centerline x"
                                              ],
                                              
                                "thickness"  : 1,
                                "separation" : "rel ly 0.02"
                    
                    }

}
            

dt = 1e-15
T = 100000*dt

#integrator = integrators.EulerCramer(dt)
integrator = integrators.VelocityVerlet(dt)

app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
app.run(T)



#TODO PARSER

#Start with fixed: Do positions, start with filling row by row and always check if occupied (brute force loop all fixed (within thickness of another) start with walls and loop iterate thickness in all other direction)

"""
Parser:
-placement from rectMesh
-placement from random
-unify without stickies split?
-work with events.
-check constraints. autofix them?
-Get working with different species.

Simulator:
-Verify solver.
-Add thermostat.
-Include measurements. Very OO.
-toFile and additional simulator options.
"""

