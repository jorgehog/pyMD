
from matplotlib.pylab import scatter, show, clf, close, ylim, xlim
import time, sys

from __init__ import controlParameters

class MDApp():
    
				def __init__(self, dt, mesh, atoms, mixTable, integrator):

  
							self.dt = dt
			        
							self.mesh = mesh
							self.ensemble = atoms
			       
							self.integrator = integrator
							self.ensemble.initialize(mixTable, self.mesh)
	        
        
        
				def run(self, T):
        

							n = int(T/self.dt)
							
         
							for i in range(n):
							
							#			self.makeDCVizFile()
							#			self.DCViz.mainloop()
											self.ensemble.calculateForces()
   
#   
#   
#											self.makeDCVizFile()
#											if i == 100:
#													sys.exit(0)
#											M = 0
#											for atom in self.ensemble.atoms:
#													m = max(abs(atom.pos[0]), abs(atom.pos[1]))
#													scatter(atom.pos[0], atom.pos[1])
#											  
#													M = max(M, m)
#																
#											print "max=", M
#											ylim([0,4])
#											xlim([0,2])
#											show()
#											time.sleep(1)
#											close()
#											clf()
					
											print "%d%%" % int((i+1.) / n * 100)        											
											self.integrator.updateAtoms(self.ensemble.atoms)
											self.ensemble.getKineticEnergy()      

											           

	
	

				def makeDCVizFile(self):

#		sx = ""
#		sy = ""
#		sz = ""
#		for atom in self.ensemble.atoms:
#			sx += "%g " % atom.pos[0]
#			sy += "%g " % atom.pos[1]
#			
#			if self.mesh.dim == 3:
#				sz += "%g " % atom.pos[2]
#		
#		s = [sx, sy]
#		if self.mesh.dim == 3:
#			s.append(sz)
#			
#		out = (("%g "*self.mesh.dim).strip() + "\n%s"*self.mesh.dim) % (tuple(self.mesh.shape)+ tuple(s))


								out = "%s %s\n" % tuple(self.mesh.shape)
          
								for atom in self.ensemble.atoms:
										out += "%g %g\n" % (atom.pos[0], atom.pos[1])

								with open("./MD_out.dat", "w") as f:
										f.write(out)
		