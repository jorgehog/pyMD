
from numpy import empty, random as npRandom, array

from random import uniform



class atom():
        
    def initialize(self, mesh, sigma, eps, mass, sticky, initPos, initVel):
        
        self.sigma, self.eps, self.mass = sigma, eps, mass
        
        self.sticky = sticky
    
        self.pos = initPos
        self.vel = initVel
        self.force = empty(mesh.dim)
        
        self.mesh = mesh       

    def incrementVel(self, dV):
        
        self.vel += dV
    
    def incrementPos(self, dR):
        
        self.pos = self.mesh.checkNewPos(self.pos + dR)
        
    def updateForce(self, newForce):
        
        self.force += newForce
        
    def resetForce(self):
        
        self.force.fill(0)

import re
class parser():
    
    anyNumber = "\-?\d+\.?\d*[eE]?[+-]?\d*\.?\d*"
    knownSubs = {
                 "lx" : "self.mesh.lx", 
                 "ly" : "self.mesh.ly", 
                 "lz" : "self.mesh.lz"
    }
    
    parsedData = {}
    parsedDataConstraints = {}
    parsedDataEvents = {}

    def __init__(self, mesh):
        
        self.mesh = mesh        
        
        self.types = ["free", "fixed"]
        
        self.randomMethods = ["uniform", "normal"]

        self.meshSpecs = ["boundary l|r|t|b|f|b", "centerline x|y|z"]  
        
        self.keys  = {
                    "fraction"  : ["number%sum=1",
                                   "automatic:calculateAfter", 
                                   "remaining:calculateDifference"],
                    "nSpecies"  : ["number"],
                    "sigmas"    : ["numbers%len=nSpecies"],
                    "epses"     : ["numbers%len=nSpecies"],
                    "masses"    : ["numbers%len=nSpecies"],
                    
                    "partition" : ["cyclic"],
                    "positions" : ["random __method__", 
                                   "list%elements=meshSpecs"],
                    "velocities": ["random __method__"],
                     
                    "thickness" : ["number%position!=random"],
                    "separation": ["number%position!=random"],
                    
        } 
    
    
    def unpackMix(self, mixProp):
        
        for key in mixProp.keys():
            
            if key not in self.types:
                raise Exception
                
            
            self.parsedData[key] = {}
            self.parsedDataConstraints[key] = {}
            self.parsedDataEvents[key] = {}
            
            pdk = self.parsedData[key]
            pdkc = self.parsedDataConstraints[key]
            pdke = self.parsedDataEvents[key]
            
            #Loops through all keys
            for parsekey in mixProp[key].keys():
                
                #Get the correct option and constraints
                pdk[parsekey], pdkc[parsekey], pdke[parsekey] = \
                                self.getInfo(parsekey, mixProp[key][parsekey])
    
                #apply relfunc and asfunc     
                print pdk[parsekey]                           
                if pdk[parsekey] == "numbers":
                    for i, number in enumerate(mixProp[key][parsekey]):
                        if type(number) is str:
                            if "as" in number:
                                mixProp[key][parsekey][i] = self.asFunction(number, mixProp[key][parsekey])
                            elif "rel" in number:
                                mixProp[key][parsekey][i] = self.relFunction(number, mixProp[key][parsekey])
                                
                else:
                    if type(mixProp[key][parsekey]) is str:
                        if "rel" in mixProp[key][parsekey]:
                            mixProp[key][parsekey] = self.relFunction(mixProp[key][parsekey])
        
        for key, val in mixProp.items():
            for parsekey, parseval in val.items():
                if type(parseval) is str and self.parsedData[key][parsekey] == "number":
                    if "as" in parseval:
                        target = parseval.lstrip("as ")
                        
                        if target in self.types:
                            targetval = mixProp[target][parsekey]
                            
                            if type(targetval) is str:
                                if "as" in targetval:
                                    raise Exception
                                
                            mixProp[key][parsekey] = targetval
                    
                elif self.parsedData[key][parsekey] == "numbers":
                    if "as" in parseval:
                        target = parseval.lstrip("as ")
                        
                        if target in self.types:
                            targetval = mixProp[target][parsekey]
                            
                            if type(targetval) is str:
                                if "as" in targetval:
                                    raise Exception
                                
                            mixProp[key][parsekey] = targetval
                
            
        for k, v in mixProp.items():
            print k, ":"
            for _k, _v in v.items():
                print _k, "  ", _v                        
                        
                        
        return
        
    def checkOpt(self, opt, value):

        #Direct match
        if type(value) is str:
                    
            if value in opt:
                return True
            
        #if number is allowed, check if value is a number
        if "number" in opt and "numbers" not in opt:
            if type(value) in [float, int]:
                return True
            elif type(value) is str:
                if "as" in value or "rel" in value:
                    return True
        
        #and so on
        if "numbers" in opt:
            if type(value) == list:
                for val in value:
                    if type(val) in [float, int]:
                        return True
            elif type(value) is str:
                if "as" in value or "rel" in value:
                    return True
        
        if "list" in opt:
            if type(value) == list:
                return True
                
        return False
        
    def getInfo(self, key, value):
        
        options = self.keys[key]
        
        option = ""
#        print "-------|---------"
#        Loop over all listed allowed options for given key
        for opt in options:
            
            if type(value) is str:
                for method in self.randomMethods:
                    if method in value:
                        opt = opt.replace("__method__", method)            
            
#            print value, opt, type(value), type(opt)
            #Check if value is in allowed options
            if self.checkOpt(opt, value):
#                print "found option ", opt
                option = opt
                break
       
        #Get constraints and events
        constraint = event = rawOption = None
        if "%" in option:
            rawOption, constraint = option.split("%")
        
        elif ":" in option:
            rawOption, event = option.split(":")
        else:
            rawOption = option
        
#        print rawOption, constraint, event
#        print "----"
        return rawOption, constraint, event
            
            
            
        
    
    def getNumbers(self, _list):
        
        return
        
    def asFunction(self, expr, _list = None):
        
        i =  int(re.findall("as\s(\d+)", expr)[0])
        
        return _list[i]
        
    def relFunction(self, expr, _list=None):
        
        
        print "---"
        print expr
        pattern = "rel\s(%s)\s(%s)" % ("%s|%s" % (self.anyNumber, "|".join(self.knownSubs.keys())), self.anyNumber)
        print pattern
        _id, scale = re.findall(pattern, expr)[0]
        
        scale = float(scale)
        
        print _id, scale
        if _id and _id in self.knownSubs.keys(): #parameter
            
            return eval(self.knownSubs[_id])*scale
       
       
        elif _id and type(eval(_id)) in [int]: #index
            
            return _list[int(_id)]*scale
        
        
            

        

class ensemble():
    
    def __init__(self, N, forceModel):
        
        self.atoms = [atom() for i in range(N)]
        self.forceModel = forceModel
        self.N = N
        
        
    def initialize(self, mixingProperties, mesh):
        
        self.mesh = mesh
        mixTable = self.unpackMix(mixingProperties)  
        
        
        for i, _atom in enumerate(self.atoms):
            _atom.initialize(mesh, *mixTable[i])
            
        
    
    def unpackMix(self, mixProp):
        
        unpacker = parser(self.mesh)
        
        unpacker.unpackMix(mixProp)
        
        import sys
        sys.exit(0)
        
        return [[1, 10, 1, False, array([uniform(0, self.mesh.lx), uniform(0, self.mesh.ly)]), 0.1*npRandom.normal(size=self.mesh.dim)] for i in range(self.N)] #TODO
        
    def getRelPos(self, atom1, atom2):
        #SMARTIFY THIS
    
        minRelPos = atom1.pos - atom2.pos
        minRelPos2 = (minRelPos**2).sum()
        

     
        for i in range(self.mesh.dim):
            if self.mesh.p[i]:

                tmp = atom2.pos[i]
                atom2.pos[i] -= self.mesh.shape[i]
                
                test = ((atom1.pos - atom2.pos)**2).sum()
                if test < minRelPos2:
                    minRelPos = atom1.pos - atom2.pos
                    minRelPos2 = test
                    
                atom2.pos[i] += 2*self.mesh.shape[i]
                
                test = ((atom1.pos - atom2.pos)**2).sum()
                if test < minRelPos2:
                    minRelPos = atom1.pos - atom2.pos
                    minRelPos2 = test
                    
                atom2.pos[i] = tmp
                
        return minRelPos, minRelPos2
        
        
    def calculateForces(self):
        
        for atom in self.atoms:
            atom.resetForce()
        
        for i, atom1 in enumerate(self.atoms[:-1]):         
            for k, atom2 in enumerate(self.atoms[(i+1):]):
                if atom1.sticky and atom2.sticky:
                    continue
                
                relPos, relPos2 = self.getRelPos(atom1, atom2)                
                
                force = self.forceModel.calculateForce(atom1, atom2, relPos, relPos2)
                
                atom1.updateForce(force)
                atom2.updateForce(-force)
        