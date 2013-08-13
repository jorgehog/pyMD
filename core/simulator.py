
import os, signal, glob, shutil
from os.path import join

from DCVizWrapper import DCVizThread

from __init__ import controlParameters

class MDApp():
    
    def __init__(self, dt, mesh, atoms, mixTable, integrator):
        
        self.dt = dt
        self.mesh = mesh
        self.ensemble = atoms
        self.integrator = integrator
  
        self.ensemble.setSimulator(self)
        self.ensemble.initialize(mixTable, self.mesh)
        
        self.cwd = os.getcwd()            
        self.cleanFiles()
        
        self.makeDCVizFile()        
        self.DCVizApp = DCVizThread(join(self.cwd, "MD_out0.dat"), True, False, delay = 0.1)
        
        
        signal.signal(signal.SIGINT, self.signal_handler)        
        
        self.stopped = False
        
        
        
    def cleanFiles(self):
        
        files = glob.glob("MD_out*.dat")

        for _file in [join(self.cwd, __file) for __file in files]:
            os.remove(_file)
        
        
    
    def run(self, T):
        
        n = int(T/self.dt)
        self.DCVizApp.start() 
        
        i = 0
        while i < n and not self.stopped:

            self.integrator.updateAtoms(self.ensemble)

            self.ensemble.getKineticEnergy()
            self.ensemble.checkLinearMomentum()
            
            if (i % 100 == 0):
                self.makeDCVizFile(i)
                print "%d%%" % int((i+1.) / n * 100)   

            i += 1
            
        self.clean()
        
        
        
    def clean(self):    
        self.DCVizApp.stop()


    def makeDCVizFile(self, i=0):
        
        out = "%s %s\n" % tuple(self.mesh.shape)    
        
        for atom in self.ensemble.atoms:
            out += "%g %g\n" % (atom.pos[0], atom.pos[1])

        with open(join(self.cwd, "MD_out%d.dat" % i), "w") as f:
            f.write(out)


    def signal_handler(self, signal, frame):

        self.stopped = True
        self.clean()
        