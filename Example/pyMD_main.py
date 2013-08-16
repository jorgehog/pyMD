
from pyMD.core import simulator, forceModels, integrators, thermostats

from pyMD.factory import atoms, geometry

from pyMD.misc.RNGFunctions import randomOnMesh, normalFromTermperature, allZero


shape = [1E-9, 2*1E-9]
periodicity = [False, True]
T0 = 119.74
N =150

mesh = geometry.rectMesh(periodicity, shape)

forceModel = forceModels.LennardJones()
#forceModel = forceModels.CoulombLike()

atoms = atoms.ensemble(N, forceModel, T0)

particleMixTable = {

                    "free" : {
                                "fraction"  : "remaining",
                                "nSpecies"  : 2,
                                "sigmas"    : [3.405E-10, "rel 0 0.99"],
                                "epses"     : [119.8*1.38E-23, "rel 0 0.99"],
                                "masses"    : [18*1E-27, "rel 0 0.99"],
                                
                                "positions" : randomOnMesh(),
#                                "positions" : ["boundary b"],
                                "velocities": normalFromTermperature(T0),
#                                "velocities": allZero(),
#                                
                                "partition" : "random",
                                    
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
                                              
                                "thickness"  : 1,
                                "separation" : "rel ly 0.02"
                    
                    }

}
            

dt = 1e-15
T = 1000*dt

#integrator = integrators.EulerCramer(dt)
integrator = integrators.VelocityVerlet(dt)

#thermostat = thermostats.BerendsenThemostat(T0, dt, 15*dt)
thermostatTop = thermostats.BerendsenThemostat(T0*1.5, dt, 15*dt)
thermostatBot = thermostats.BerendsenThemostat(T0*1.5, dt, 15*dt)
thermostatMid = thermostats.BerendsenThemostat(T0*0.5, dt, 15*dt)

app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
#app.addField(thermostat, geometry.predefinedRegions("allMesh", mesh=mesh))
app.addField(thermostatTop, geometry.predefinedRegions("top", mesh=mesh, width=shape[1]/5.0))
app.addField(thermostatBot, geometry.predefinedRegions("bottom", mesh=mesh, width=shape[1]/5.0))
app.addField(thermostatMid, geometry.predefinedRegions("midX", mesh=mesh, width=shape[1]/5.0))

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

