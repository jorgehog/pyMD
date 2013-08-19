
import re

from math import sqrt
from numpy import array, zeros, string_
from random import randint        

from pyMD.misc import RNGFunctions
        
        

class plcDist():
    
    nav = {
           "l" : {
               "init" : array([0, 0]),
               "iter" : array([0, 1]), 
               "next" : array([1, 0]),
               "end"  : array([0, "ly"])               
           },
           
           "r" : {
               "init" : array(["lx", 0]), 
               "iter" : array([0, 1]), 
               "next" : array([-1, 0]),
               "end"  : array(["lx", "ly"])
           },
    
           "t" : {
               "init" : array([0, "ly"]), 
               "iter" : array([1, 0]), 
               "next" : array([0, -1]),
               "end"  : array(["lx", "ly"])
           },
           
           "b" : {
               "init" : array([0, 0]), 
               "iter" : array([1, 0]), 
               "next" : array([0, 1]),
               "end"  : array(["lx", 0])
           },
           
           "x" : {
               "init" : array([0, "ly/2. - width/2."]), 
               "iter" : array([1, 0]), 
               "next" : array([0, 1]),
               "end"  : array(["lx", "ly/2. - width/2."])
           },
           
           "y" : {
               "init" : array(["lx/2. - width/2.", 0]), 
               "iter" : array([0, 1]), 
               "next" : array([1, 0]),
               "end"  : array(["lx/2. - width/2.", "ly"])
           }
               
    }    
    
    allPos = []
    
    def __init__(self, mesh, thickness, separation, *positions):
        
        self.mesh = mesh
        self.thickness = thickness
        self.separation = separation
        self.finished = False
        
        #Make a copy of all the keys and initialize empty lists
        self.positions={}
        self.allKeys = tableParser.meshSpecs.keys()        
 
        for pos in positions:
            place, align = pos.split()
            
            if align not in tableParser.meshSpecs[place]:
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
                        for __key in tableParser.knownSubs.keys():   
              
                            if __key in element:
                                element = element.replace(__key, tableParser.knownSubs[__key])
                
                    if doEval:
                        tmp[i] = eval(element)
                    else:
                        tmp[i] = element

                self.nav[key][_key] = tmp.copy()                        
                        

            
            
        self.curStep = 0
        self.keyIndex = 0
        self.subKeyIndex = 0
        
        
    def isOccupied(self, pos):
        
        posWithBounds = self.mesh.checkNewPos(pos)
        thresh = self.separation*(1 - 1E-3)       
        
        for setPos in self.allPos:
            
            if self.arrayLength(posWithBounds - setPos) < thresh or self.arrayLength(pos - setPos) < thresh:
#                print "Occupied at pos: ", pos, self.arrayLength(pos - setPos), "from", setPos
                return True
                
        return False
        
    
    def arrayLength(self, a):
        
        return sqrt((a**2).sum())
    
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
        
        
        
        key = self.positions[self.keys[self.keyIndex]][self.subKeyIndex]
        nav = self.nav[key]
        
        #Calculate the lineSize and width in indexes
        lineSize = int(self.arrayLength(nav["end"] - nav["init"])/self.separation) + 1 #Number of indexes on a line
        
        #Curindex = widthPos*lineSize + linePos
        linePos = self.curStep%lineSize #Position along the line
        widthPos = (self.curStep - linePos)/lineSize #Position along the width
        
        lastIndex = self.thickness*lineSize - 1
        
#        print linePos, "out of", lineSize        
#        print widthPos, "out of", self.thickness       
#        print self.curStep, "out of", lastIndex
        
        #Iterate the indexes and keys
        if self.curStep == lastIndex:
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
        
        if self.isOccupied(pos):
            return None
        else:
            self.allPos.append(pos)          
            return pos
        
        
            
            

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
        except:
            thickness = 1
            
        return positions, velIter, separation, thickness
    
    
    
    def getPosAndVelIterators(self, props):
        
    
        positions, velIter, separation, thickness = self.getVelAndPosSpecsFromProps(props)
                
        
        if type(positions) != list:
            if type(positions) is str:
                
                posIter = plcDist(self.mesh, thickness, separation, positions)
                
            elif isinstance(positions, RNGFunctions.RNGs):
                
                posIter = positions
                    
            else:
    
                raise Exception("Undefined position argument", positions)                
                
                
        else:
            posIter = plcDist(self.mesh, thickness, separation, *positions)
    
    
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
                
        n = 0
        
        while n < nCap:
            
            sigma, eps, mass = self.getObjectProps(props, partFunc, n)            

            print nCap, n, self.nInitialized, self.N
            
            newPos = None
            while newPos is None and not posIter.finished:
                newPos = posIter(sigma=sigma, eps=eps, mass=mass, ensemble=self.ensemble, dim=self.mesh.dim)
            
            if posIter.finished:
                break
            
            newVel = velIter(sigma=sigma, eps=eps, mass=mass, ensemble=self.ensemble, dim=self.mesh.dim)
            
            
            self.ensemble.atoms[self.nInitialized].initialize(self.mesh, 
                                                              self.N,
                                                              sticky, 
                                                              sigma, 
                                                              eps, 
                                                              mass, 
                                                              newPos, 
                                                              newVel)
                                                              
#            table.append([sticky, sigma, eps, mass, newPos, newVel])
            
#            print "n=", n
            self.nInitialized += 1
            n += 1
#        print "done"
            
        
#        if appendTo:
#            return appendTo + table
#        return table
        
        
        
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