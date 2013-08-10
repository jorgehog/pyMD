

class rectMesh():
    
    lz = None
    pz = None
    
    def __init__(self, periodicity, shape):
        
        self.px, self.py = shape[:2]     
        self.dim = len(shape)
        
        if self.dim == 3: 
            self.pz = shape[2]
            
        self.p = [self.px, self.py, self.pz]
        
        self.setShape(shape)


    def setShape(self, shape):
   
        self.shape = shape
   
        self.lx, self.ly = shape[:2]        
        
        if self.dim == 3:
            self.lz = shape[2]
            
            
    def checkNewPos(self, pos):

        for i in range(self.dim):
            
            if pos[i] > self.shape[i]:
                if self.p[i]:
                    pos[i] = pos[i]%self.shape[i]
                else:
                    pos[i] = self.shape[i]
                    #Reverse transverse momentum?
                    
            elif pos[i] < 0:
                if self.p[i]:
                    pos[i] = pos[i]%self.shape[i]
                else:
                    pos[i] = 0
                    #Reverse transverse momentum?

        return pos
        