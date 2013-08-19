

class rectMesh:
    
    lz = None
    pz = None
    
   
    normals3D = [] #...    
    
    def __init__(self, periodicity, shape):
        
        self.px, self.py = periodicity[:2]     
        self.dim = len(shape)
        
        if self.dim == 3: 
            self.pz = periodicity[2]
            self.normals = self.normals3D
        else:
            self.normals = lambda i: i==0 and 1 or 0
            
        self.p = [self.px, self.py, self.pz]
        
        self.setShape(shape)

    def normals3d(self):
        pass

    def setShape(self, shape):
   
        self.shape = shape
        self.volume = shape[0]*shape[1]        
        
        self.lx, self.ly = shape[:2]        
        
        if self.dim == 3:
            self.lz = shape[2]
            self.volume*=shape[2]
            
            
    def checkNewPos(self, pos, vel=None):

        for i in range(self.dim):
                      
            if pos[i] >= self.shape[i]:
                if self.p[i]:
                    pos[i] = pos[i]%self.shape[i]
                else:
                    if vel is not None:
                        vel[i] *= -1
                    #Reverse transverse momentum?
                    
            elif pos[i] < 0:
                if self.p[i]:
                    pos[i] = pos[i]%self.shape[i]
                else:
                    if vel is not None:
                        vel[i] *= -1
        
        return pos



class region:

    topology = []    
    atoms = []
    
    def __init__(self, x0, x1=None, y0=None, y1=None, z0=None, z1=None, description = None):
        
        if x1 is not None and (y0 is None or y1 is None):
            raise Exception()
        
        if type(x0) is list:
            
            if x1 or y0 or y1 or z0 or z1:
                raise Exception()
            
            self.topology = x0
        else:
            if x1 is None:
                raise Exception()
    
            self.topology = [[x0, x1], [y0, y1]]
        
        if z0 is not None and z1 is not None:
            self.topology.append([z0, z1])
    
        self.dim = len(self.topology)
        
        self.volume = 1
        for axis in self.topology:
            self.volume *= axis

        if not description:
            description = "V = %g | loc = %s" % (self.volume, str(tuple(self.topology)))
            
        self.description = description

    def __eq__(self, other):
        
        return self.topology == other.topology
        
    def __str__(self):
        
        return "<region> %s" % self.description

    def append(self, atom):
                
        for i in range(self.dim):
            if atom.pos[i] < self.topology[i][0] or atom.pos[i] > self.topology[i][1]:
                return


        self.atoms.append(atom)                
                
    def reset(self):
#        print id(self.atoms)
        self.atoms = []
#        print id(self.atoms)
#        print "----"
        
    

def predefinedRegions(name, **kwargs):
    
        allMesh = [
                    [0, 'kwargs["mesh"].lx'],
                    [0, 'kwargs["mesh"].ly']
            ]
            
        top = [
                [0, 'kwargs["mesh"].lx'],
                ['kwargs["mesh"].ly - kwargs["width"]', 'kwargs["mesh"].ly']        
        
        ]
        
        bottom = [
                [0, 'kwargs["mesh"].lx'],
                [0, 'kwargs["mesh"].ly + kwargs["width"]']      
        
        ]
    
        midX = [
                [0, 'kwargs["mesh"].lx'],
                ['kwargs["mesh"].ly/2.0 - kwargs["width"]/2.0', 'kwargs["mesh"].ly/2.0 + kwargs["width"]/2.0']        
        
        
        ]
        
        allPredefines = [allMesh, top, bottom, midX]
         
        for s, preDef in enumerate(allPredefines):
            for i, _l in enumerate(preDef):
                print _l
                for j, _k in enumerate(_l):
                    allPredefines[s][i][j] = eval(str(_k))
                    
        return region(eval(name), description=name)
        
                
        
        

        

        
        
        
    
    
    


