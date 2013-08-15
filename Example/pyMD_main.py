
from pyMD.core import simulator, forceModels, integrators

from pyMD.factory import atoms, geometry

from pyMD.misc.RNGFunctions import randomOnMesh, normalFromTermperature, allZero


shape = [1E-9, 2*1E-9]
periodicity = [False, True]
T = 100
N = 250

mesh = geometry.rectMesh(periodicity, shape)

forceModel = forceModels.LennardJones()
#forceModel = forceModels.CoulombLike()

atoms = atoms.ensemble(N, forceModel)

particleMixTable = {

                    "free" : {
                                "fraction"  : "remaining",
                                "nSpecies"  : 2,
                                "sigmas"    : [3.405E-10, "rel 0 0.99"],
                                "epses"     : [119.8*1.38E-23, "rel 0 0.99"],
                                "masses"    : [18*1E-27, "rel 0 0.99"],
                                
                                "positions" : randomOnMesh(),
                                "velocities": normalFromTermperature(T),
#                                "velocities": allZero(),
                                
                                "partition" : "cyclic"
                                    
                    },
                             
                             
                    "fixed" : {
    
                                "fraction"  : "automatic",
#                                "fraction" : 0,
                                
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
                                              
                                "thickness"  : 2,
                                "separation" : "rel ly 0.02"
                    
                    }

}
            

dt = 1e-16
T = 1000*dt

#integrator = integrators.EulerCramer(dt)
integrator = integrators.VelocityVerlet(dt)

app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
app.run(T)



#TODO PARSER

#Start with fixed: Do positions, start with filling row by row and always check if occupied (brute force loop all fixed (within thickness of another) start with walls and loop iterate thickness in all other direction)

"""
Parser:
-work with events.
-check constraints. autofix them?
-Get working with different species.

Simulator:
-Verify solver.
-Add thermostat.
-Include measurements. Very OO.
-toFile and additional simulator options.
"""

