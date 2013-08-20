
import re

from math import pi
from random import randint        

from pyMD.misc import RNGFunctions
        
        


        
        
            
            

class tableParser:
    
    anyNumber = "\-?\d+\.?\d*[eE]?[+-]?\d*\.?\d*"
    knownSubs = {
                 "lx" : "self.mesh.lx", 
                 "ly" : "self.mesh.ly", 
                 "lz" : "self.mesh.lz",
                 "width" : "self.thickness*self.separation"
    }
    
    
    types = ["free", "fixed"]
        
    randomMethods = ["uniform", "normal"]

    meshSpecs = {"boundary" : ["l","r","t","b","f","b"], "centerline" : ["x", "y", "z"]}  
        
    keys  = {
                "fraction"  : ["number%sum=1",
                               "automatic:calculateAfter", 
                               "remaining:calculateDifference"],
                "nSpecies"  : ["number"],
                "sigmas"    : ["numbers%len=nSpecies"],
                "epses"     : ["numbers%len=nSpecies"],
                "masses"    : ["numbers%len=nSpecies"],
                
                "partition" : ["cyclic", "random"],
                "positions" : ["random __method__", 
                               "list%elements=meshSpecs"],
                "velocities": ["random __method__"],
                 
                "thickness" : ["number%position!=random"],
                "separation": ["number%position!=random"],
                
    }     
    
    
    parsedData = {}
    parsedDataConstraints = {}
    parsedDataEvents = {}
    
    nInputAtoms = 6

    def __init__(self, ensemble, mesh, N):
        
        self.mesh = mesh        
        
        self.ensemble = ensemble
        
        self.N = N        
        
        self.nInitialized = 0
        
        self.volumeUsed = 0
        
        self.occupied = {}
        
 
    def setEnsembleValues(self, prop):

        freeSigmas = fixedSigmas = None
        
        if "fixed" in prop.keys():
            fixedSigmas = prop["fixed"]["sigmas"]

        if "free" in prop.keys():
            freeSigmas = prop["free"]["sigmas"]
            
        setattr(self.ensemble, "fixedSigmas", fixedSigmas)
        setattr(self.ensemble, "freeSigmas", freeSigmas)
        
        
                    
            
    
    def unpackMix(self, mixProp):
        
        mixProp = self.getValuesConstraitsEvents(mixProp)
        
        self.setEnsembleValues(mixProp)
        
        self.initializeTable(mixProp)
        
        
    def initializeTable(self, mixProp):
        
        self.unpackProps(mixProp["fixed"], sticky = True)
        setattr(self.ensemble, "nFixed", self.nInitialized)
        
        self.unpackProps(mixProp["free"], sticky = False)
        setattr(self.ensemble, "nFree", self.nInitialized - self.ensemble.nFixed)

    

    def getNCapFromProps(self, props):
        
            fraction  = props["fraction"]    
              
            if fraction == "automatic" or fraction == "remaining":
                m = self.N - self.nInitialized
                
                return (m > 0 and m) or 0

                    
            elif type(fraction) in [float, int]:
                    
                if 1 < fraction or fraction < 0: 
                    raise Exception("Fraction must be on [0, 1], got %g" % fraction)
                
                return int(self.N*fraction)
                
            else: 
                raise Exception("Unknown fraction model: ", fraction)
    
    
    
    def getVelAndPosSpecsFromProps(self, props):
        
        try:
            positions = props['positions']
        except:
            positions = RNGFunctions.randomOnMesh(self.mesh)
        
        try:
            velIter = props['velocities']
        except:
            velIter = RNGFunctions.allZero()
        
        try:
            separation = props['separation']
        except:
            separation = min(self.mesh.shape)/self.N
        
        try:
            thickness = props['thickness']
            
            if thickness == "full":
                thickness = max(self.mesh.shape)/separation
                
            
        except:
            thickness = 1
            
        return positions, velIter, separation, thickness
    
    
    
    def getPosAndVelIterators(self, props):
        
    
        positions, velIter, separation, thickness = self.getVelAndPosSpecsFromProps(props)
                
        
        if type(positions) != list:
            if type(positions) is str:
                
                self.occupied[positions] = thickness*separation                            
                
                posIter = RNGFunctions.plcDist(self.mesh, thickness, separation, self.meshSpecs, self.knownSubs, 0, positions)
                
            elif isinstance(positions, RNGFunctions.RNGs):
                
                posIter = positions
                    
            else:
    
                raise Exception("Undefined position argument", positions)                
                
                
        else:
            
            for pos in positions:
                    self.occupied[pos] = thickness*separation            
            
            posIter = RNGFunctions.plcDist(self.mesh, thickness, separation, self.meshSpecs, self.knownSubs,  0,  *positions)
    
    
        setattr(posIter, "_separation", separation)
        
        return posIter, velIter
    
    def getPartFuncFromProps(self, props):
        
        try:
            nSpecies = props["nSpecies"]
        except:
            nSpecies = 1
        
        try:
            partition = props["partition"]
        except:
            partition = "random"
        
        if partition == "cyclic":
            partFunc = lambda i: i%nSpecies
        elif partition == "random":
            partFunc = lambda i: randint(0, nSpecies-1)
        else:
            raise Exception("Unknown partition function", partition)
        
        return partFunc
            
    
    def getObjectProps(self, props, partFunc, n):
        
        i = partFunc(n)
        
        sigma = props["sigmas"][i]
        eps = props["epses"][i]
        mass = props["masses"][i] 
    
        return sigma, eps, mass    
    
    def unpackProps(self, props, sticky):
    
        nCap = self.getNCapFromProps(props)
        
        if nCap == 0:
            return
        
        posIter, velIter = self.getPosAndVelIterators(props)
        partFunc = self.getPartFuncFromProps(props)
        
        posIter.prep(N=nCap, parser=self)
        velIter.prep(N=nCap, parser=self)
                
        n = 0
        
        while n < nCap:
            
            sigma, eps, mass = self.getObjectProps(props, partFunc, n)            

            print nCap, n, self.nInitialized, self.N
            
            newPos = None
            while newPos is None and not posIter.finished:
                newPos = posIter(sigma=sigma, eps=eps, mass=mass, parser=self)
#                print posIter.finished

            
            if posIter.finished:
                print "unable to place more than %d atoms" % n
                break
            
            newVel = velIter(sigma=sigma, eps=eps, mass=mass, parser=self)
            
            
            self.ensemble.atoms[self.nInitialized].initialize(self.mesh, 
                                                              self.N,
                                                              sticky, 
                                                              sigma, 
                                                              eps, 
                                                              mass, 
                                                              newPos, 
                                                              newVel)
                                                              

            self.nInitialized += 1
            n += 1
        
        

        
        
        
    def getValuesConstraitsEvents(self, mixProp):
        for key in mixProp.keys():
            
            if key not in self.types:
                raise Exception("%s is not an allowed type." % key)
                
            
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
#                print pdk[parsekey]                           
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
                                    raise Exception("Loop in 'as' statement.")
                                
                            mixProp[key][parsekey] = targetval
                    
                elif self.parsedData[key][parsekey] == "numbers":
                    if "as" in parseval:
                        target = parseval.lstrip("as ")
                        
                        if target in self.types:
                            targetval = mixProp[target][parsekey]
                            
                            if type(targetval) is str:
                                if "as" in targetval:
                                    raise Exception("Loop in 'as' statement.")
                                
                            mixProp[key][parsekey] = targetval     
                            
        return mixProp
        
        
        
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
        
        
#        print "---"
#        print expr
        pattern = "rel\s(%s)\s(%s)" % ("%s|%s" % (self.anyNumber, "|".join(self.knownSubs.keys())), self.anyNumber)
#        print pattern
        _id, scale = re.findall(pattern, expr)[0]
        
        scale = float(scale)
        
#        print _id, scale
        if _id and _id in self.knownSubs.keys(): #parameter
            
            return eval(self.knownSubs[_id])*scale
       
       
        elif _id and type(eval(_id)) in [int]: #index
            
            return _list[int(_id)]*scale