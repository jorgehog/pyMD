
import sys, os#, threading

cwd = os.getcwd()

sys.path.append(cwd + "/../../")
#sys.path.append("E:/programmer/DCViz-master/src")
#
#with open("MD_out.dat", "w") as f:
#    f.write("")
#
#from DCViz_classes import MD_out
#
#class startDCVIZ(threading.Thread):
#    
#    MD = MD_out(os.path.join(cwd, "MD_out.dat"), True)
#    
#    def run(self):
#        
#        self.MD.mainloop()
#        
#    def stop(self):
#        self.MD.stopped = True
#        print "stopped viz thead"
        
        
        
        

from pyMD.core import simulator, forceModels, integrators

from pyMD.factory import atoms, geometry


shape = [2, 4]
periodicity = [False, True]
N = 100

mesh = geometry.rectMesh(periodicity, shape)

forceModel = forceModels.LennardJones()
#forceModel = forceModels.CoulombLike()

atoms = atoms.ensemble(N, forceModel)

particleMixTable = {

                    "free" : {
                                "fraction"  : "remaining",
                                "nSpecies"  : 2,
                                "sigmas"    : [0.01, "rel 0 1e-5"],
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
            

dt = 1e-15
T = 10*dt

integrator = integrators.EulerCramer(dt)

#visualizer = startDCVIZ()
#visualizer.start()

app = simulator.MDApp(dt, mesh, atoms, particleMixTable, integrator)
app.run(T)

#visualizer.stop()







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

