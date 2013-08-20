

from math import sqrt, pi
from numpy import array, zeros, string_
from numpy.random import uniform, normal
from scipy.constants import k as kB


class RNGs:
    
    finished = False
    
    def prep(self, **kwargs):
        return
    
    def isOccupied(self, ens, pos, thresh):
        
        a = RNGs()
        setattr(a, "pos", pos)
        
        for atom in ens.atoms:
            if atom.initialized:
                rel, rel2 = ens.getRelPos(atom, a)
                
                if sqrt(rel2) < thresh:
#                    print pos, "too close to ", atom.pos, "distance", sqrt(rel2)
#                    print "sep = ", thresh
                    return True
        
        return False
        
    def arrayLength(self, a):
        
        return sqrt((a**2).sum())
        

class randomOnMesh(RNGs):
    
        def getNewPos(self, shape):
            
             plc = []
             
             for l in shape:
                 plc.append(uniform(0, l))
             
             return array(plc)          


        def __call__(self, **kwargs):
            
             ens = kwargs["parser"].ensemble
             
             pos = self.getNewPos(ens.mesh.shape)

             scale = 1

             k = 0
             while self.isOccupied(ens, pos, kwargs["sigma"]*scale):             
                 pos = self.getNewPos(ens.mesh.shape)
                 
                 if k > 10:
                     scale *= 0.75
                     k = 0
                 else:
                     k += 1
             
             return pos
             

class uniformOnMesh(RNGs):
    
    #finds open areas and fills them uniformly
    def prep(self, **kwargs):
        
        self.nCap = kwargs["N"]
        self.mesh = kwargs["parser"].mesh
    
        dim = self.mesh.dim
        volume = self.mesh.volume
        remainingVolume = (volume - kwargs["parser"].volumeUsed)
#        remainingVolume = 0.825*volume
#        print remainingVolume, volume, remainingVolume/volume*100, "%"  
#        raw_input() 
        
        spacePerAtom = remainingVolume/self.nCap        

        sep = (3*spacePerAtom/(pi*(dim + 1)))**(1./dim)

        thickness = int(self.mesh.lx/sep)

        pos = "boundary l"  
  
        try:
            offset = kwargs["parser"].occupied[pos]
        except:
            offset = 0
      
        
        self.wrappedIter = plcDist(self.mesh,thickness, sep, kwargs["parser"].meshSpecs, kwargs["parser"].knownSubs, offset, pos)
        
        self.wrappedIter.prep(**kwargs)
        
        
    def __call__(self, **kwargs):
        
        return self.wrappedIter(**kwargs)
    
    
    


class normalFromTermperature(RNGs):
    
    def __init__(self, T):

        self.T = T
        self.sqrtK = sqrt(kB)
            
    def __call__(self, **kwargs):
        
        return normal(scale = self.sqrtK*sqrt(self.T/kwargs["mass"]), size = kwargs["parser"].mesh.dim)
        
        
class allZero(RNGs):
        
    def __call__(self, **kwargs):
        
        return zeros(kwargs["parser"].mesh.dim)
        


class plcDist(RNGs):
     
    
#    allPos = []
    
    def __init__(self, mesh, thickness, separation, meshSpecs, knownSubs, offset, *positions):

        self.nav = {
                       "l" : {
                           "init" : array(["offset", 0]),
                           "iter" : array([0, 1]), 
                           "next" : array([1, 0]),
                           "end"  : array(["offset", "ly"])               
                       },
                       
                       "r" : {
                           "init" : array(["lx - offset", 0]), 
                           "iter" : array([0, 1]), 
                           "next" : array([-1, 0]),
                           "end"  : array(["lx - offset", "ly"])
                       },
                
                       "t" : {
                           "init" : array([0, "ly - offset"]), 
                           "iter" : array([1, 0]), 
                           "next" : array([0, -1]),
                           "end"  : array(["lx", "ly - offset"])
                       },
                       
                       "b" : {
                           "init" : array([0, "offset"]), 
                           "iter" : array([1, 0]), 
                           "next" : array([0, 1]),
                           "end"  : array(["lx", "offset"])
                       },
                       
                       "x" : {
                           "init" : array([0, "ly/2. - width/2. + offset"]), 
                           "iter" : array([1, 0]), 
                           "next" : array([0, 1]),
                           "end"  : array(["lx", "ly/2. - width/2. + offset"])
                       },
                       
                       "y" : {
                           "init" : array(["lx/2. - width/2. + offset", 0]), 
                           "iter" : array([0, 1]), 
                           "next" : array([1, 0]),
                           "end"  : array(["lx/2. - width/2. + offset", "ly"])
                       }
                           
        }   


        
        self.mesh = mesh
        self.thickness = thickness
        self.separation = separation
        self.finished = False
        
        #Make a copy of all the keys and initialize empty lists
        self.positions={}
        self.allKeys = meshSpecs.keys()        
 
        for pos in positions:
            place, align = pos.split()
            
            if align not in meshSpecs[place]:
                raise Exception("Unable to perform placement at %s" % pos)
            else:
                if place in self.positions.keys():
                    self.positions[place].append(align)
                else:   
                    self.positions[place] = [align]
                    
        self.keys = self.positions.keys()
        
        for key, value in self.nav.items():

            for _key, _value in value.items():
                tmp = zeros(shape=_value.shape)
                
                for i, element in enumerate(_value):
                    
                    doEval = False
   
                    if type(element) is string_:
                        doEval = True
                        for __key in knownSubs.keys():   
              
                            if __key in element:
                                element = element.replace(__key, knownSubs[__key])
                
                    if doEval:
                        tmp[i] = eval(element)
                    else:
                        tmp[i] = element

                self.nav[key][_key] = tmp.copy()                        
                        

            
            
        self.curStep = 0
        self.keyIndex = 0
        self.subKeyIndex = 0
    
    def prep(self, **kwargs):
        self.ens = kwargs["parser"].ensemble
        return
        
#    def isOccupied(self, pos):
#        
#        posWithBounds = self.mesh.checkNewPos(pos)
#        thresh = self.separation*(1 - 1E-3)       
#        
#        for setPos in self.allPos:
#            
#            if self.arrayLength(posWithBounds - setPos) < thresh or self.arrayLength(pos - setPos) < thresh:
##                print "Occupied at pos: ", pos, self.arrayLength(pos - setPos), "from", setPos
#                return True
#                
#        return False
        
    

    
    def __call__(self, **kwargs):
        
        
        ### DEBUG
#        out = "%s %s\n" % tuple(self.mesh.shape)    
#        
#        for pos in self.allPos:
#            out += "%g %g\n" % (pos[0], pos[1])
#
#        with open(join(getcwd(), "MD_out1337.dat"), "w") as f:
#            f.write(out)        
#        ####
#        
#        print "----"
##        if self.keyIndex == 1:
#        raw_input()
        
        
#        print self.keyIndex, self.subKeyIndex, self.finished   
        key = self.positions[self.keys[self.keyIndex]][self.subKeyIndex]
        nav = self.nav[key]
        
        #Calculate the lineSize and width in indexes
        lineSize = int(self.arrayLength(nav["end"] - nav["init"])/self.separation) + 1 #Number of indexes on a line
        
#        Curindex = widthPos*lineSize + linePos
        linePos = self.curStep%lineSize #Position along the line
        widthPos = (self.curStep - linePos)/lineSize #Position along the width
        
        lastIndex = self.thickness*lineSize - 1
        
#        print linePos, "out of", lineSize        
#        print widthPos, "out of", self.thickness       
#        print self.curStep, "out of", lastIndex
        
        #Iterate the indexes and keys
        if self.curStep == lastIndex:

            kwargs["parser"].volumeUsed += (0.5)**self.mesh.dim*((self.thickness - 1)*lineSize + linePos)*(self.separation)**2
            self.curStep = 0
            
            if self.subKeyIndex == len(self.positions[self.keys[self.keyIndex]]) - 1:
                self.subKeyIndex = 0
                self.keyIndex += 1
                
                if self.keyIndex == len(self.keys):
                    self.finished = True
                    
            else:
                self.subKeyIndex += 1
        else:
            self.curStep += 1
                
        
        pos = nav["init"] + self.separation*(nav["iter"]*linePos + nav["next"]*widthPos)
        
#        print "trying to place", pos
        
        thresh = self.separation - 1E-3
        if self.isOccupied(self.ens, pos, thresh):
            return None
        else:
#            self.allPos.append(pos)          
            return pos