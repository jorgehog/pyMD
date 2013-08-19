
from pyMD.core import simulator, forceModels, integrators, fields

from pyMD.factory import atoms, geometry

from pyMD.misc.RNGFunctions import randomOnMesh, normalFromTermperature, allZero


shape = [1E-9, 2*1E-9]
periodicity = [False, True]
T0 = 290
N =500

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
                                
#                                "positions" : randomOnMesh(),
                                "positions" : ["boundary l"],
                                "velocities": normalFromTermperature(T0),
#                                "velocities": allZero(),
#                                
                                "partition" : "random",
                                "separation": "rel ly 0.01"
                                    
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
            

dt = 1e-15
T = 1000*dt

#integrator = integrators.EulerCramer(dt)
integrator = integrators.VelocityVerlet(dt)

#thermostat = fields.BerendsenThemostat(T0, dt, 15*dt)


app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
#app.addField(thermostat, geometry.predefinedRegions("allMesh", mesh=mesh))

top = geometry.predefinedRegions("top", mesh=mesh, width=shape[1]/5.0)
bot = geometry.predefinedRegions("bottom", mesh=mesh, width=shape[1]/5.0)
mid = geometry.predefinedRegions("midX", mesh=mesh, width=shape[1]/5.0)

thermostatTop = fields.BerendsenThemostat(T0*1.5, dt, 15*dt)
thermostatBot = fields.BerendsenThemostat(T0*1.5, dt, 15*dt)
thermostatMid = fields.BerendsenThemostat(T0*0.1, dt, 15*dt)

app.addField(thermostatTop, top)
app.addField(thermostatBot, bot)
app.addField(thermostatMid, mid)

app.addField(fields.pressureField(), top)
app.addField(fields.pressureField(), bot)
app.addField(fields.pressureField(), mid)

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

