

class rectMesh():
    
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
   
        self.lx, self.ly = shape[:2]        
        
        if self.dim == 3:
            self.lz = shape[2]
            
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
