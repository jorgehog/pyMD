
import os, signal, glob
from os.path import join
from pyMD.core import thermostats


useDCViz = False
try:
    from DCVizWrapper import DCVizThread
    useDCViz = True
except:
    pass

useDCViz = False

class MDApp():
    
    regions = []
    fields = []
    
    def __init__(self, dt, mesh, atoms, mixTable, integrator, thermostat = None):
        
        self.dt = dt
        self.mesh = mesh
        self.ensemble = atoms
        self.integrator = integrator
        
        if thermostat is None:
            thermostat = thermostats.noThermostat()
            
        self.thermostat = thermostat
  
        self.ensemble.setSimulator(self)
        self.ensemble.initialize(mixTable, self.mesh)
        
        self.cwd = os.getcwd()            
        self.cleanFiles()
        
        self.makeDCVizFile()    
        
        if useDCViz:
            self.DCVizApp = DCVizThread(join(self.cwd, "DCVizFiles", "MD_out0.dat"), True, False, delay = 0.1)
        
        
        
        signal.signal(signal.SIGINT, self.signal_handler)        
        
        self.stopped = False
        
    def addField(self, field, region):
        
        #Make sure not to perform more than necessary tests if several fields share a region
        if region not in self.regions:
            self.regions.append(region)
        else:
            region = self.regions[self.regions.index(region)]
            
        field.addRegion(region)
            
        self.fields.append(field)
        
    def cleanFiles(self):
        
        os.chdir(join(self.cwd, "DCVizFiles"))
        files = glob.glob("MD_out*.dat")

        for _file in [join(self.cwd, "DCVizFiles", __file) for __file in files]:
            os.remove(_file)
            
        os.chdir(self.cwd)
        
    def updateFields(self):

        for region in self.regions:
            region.reset()
            
        for atom in self.ensemble.atoms:
            for region in self.regions:
                region.append(atom)
        
        for field in self.fields:
            field.update()
            
            print field
    
    
    def run(self, T):
        
        n = int(T/self.dt)

        self.thermalize()        
        
        if useDCViz:
            self.DCVizApp.start() 
        
        i = 0
        
        while i < n and not self.stopped:

            self.integrator.updateAtoms(self.ensemble)
            
#            self.ensemble.getKineticEnergy()
#            self.ensemble.checkLinearMomentum()
#            self.ensemble.getTemperature()
#            self.ensemble.getPressure()

            self.updateFields()
            
        
            self.makeDCVizFile(i)
            print "%d%%" % int((i+1.) / n * 100)   

            i += 1
            
        self.clean()
        
        
        
        
    def thermalize(self):

        Tnew = self.ensemble.thermalizer.bathT
        T = 0
        
        while abs(T - Tnew)/Tnew > 1E-4:
            self.integrator.updateAtoms(self.ensemble)            
            
            T = Tnew
            Tnew = self.ensemble.thermalize()
            print Tnew, self.ensemble.thermalizer.bathT
        
        
    def clean(self): 
        if useDCViz:
                self.DCVizApp.stop()


    def makeDCVizFile(self, i=0):
        
        out = "%g %g %g\n%s\n" % (tuple(self.mesh.shape) + (self.ensemble.T,) + (str([obj.T for obj in self.fields if obj._type == "thermostat"]),))    
        
        for atom in self.ensemble.atoms:
            if atom.sticky:
                out += "%d %g %g %g %g\n" % (self.ensemble.getColor(atom), atom.pos[0], atom.pos[1], 0, 0)
            else:
                out += "%d %g %g %g %g\n" % (self.ensemble.getColor(atom), atom.pos[0], atom.pos[1], atom.vel[0], atom.vel[1])

        with open(join(self.cwd, "DCVizFiles", "MD_out%d.dat" % i), "w") as f:
            f.write(out)


    def signal_handler(self, signal, frame):

        self.stopped = True
        self.clean()
        