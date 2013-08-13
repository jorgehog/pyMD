
import re
from os.path import join
from os import getcwd
from math import sqrt
from numpy import array, float32, zeros
from numpy.random import uniform, normal

class ranDist():
    
    def __init__(self, function, *_input, **kwargs):
        
        self.function = function
        self.input = _input
        self.kwargs
        
    def __iter__(self):
        yield self.function(*self.input, **self.kwargs)
        
        
        

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
        
        for key in self.allKeys:
            self.positions[key] = []
        
        self.keys = self.positions.keys()
        
        for pos in positions:
            place, align = pos.split()
            
            if align not in tableParser.meshSpecs[place]:
                raise Exception("Unable to perform placement at %s" % pos)
            else:
                self.positions[place].append(align)
            
        
        for key, value in self.nav.items():
            print key
            for _key, _value in value.items():
                tmp = zeros(shape=_value.shape)
                
                for i, element in enumerate(_value):                    
                    if element in tableParser.knownSubs.keys():
                        
                        tmp[i] = eval(tableParser.knownSubs[element])

                    else:
                        
                        tmp[i] = element

                self.nav[key][_key] = tmp.copy()                        
                        
                print _key, self.nav[key][_key]
            
            
        self.curStep = 0
        self.keyIndex = 0
        self.subKeyIndex = 0
    
    def isOccupied(self, pos):
        
        for setPos in self.allPos:
  
            if self.arrayLength(pos - setPos) < self.separation:
                print "Occupied at pos: ", pos
                True
                
        return False
        
    
    def arrayLength(self, a):
        
        return sqrt((a**2).sum())
    
    def next(self):
        
        ### DEBUG
        out = "%s %s\n" % tuple(self.mesh.shape)    
        
        for pos in self.allPos:
            out += "%g %g\n" % (pos[0], pos[1])

        with open(join(getcwd(), "MD_out1337.dat"), "w") as f:
            f.write(out)        
        ####
        
        print "----"
        raw_input()
        
        
        
        key = self.positions[self.keys[self.keyIndex]][self.subKeyIndex]
        nav = self.nav[key]
        
        #Calculate the lineSize and width in indexes
        lineSize = self.arrayLength(nav["end"] - nav["init"])/self.separation #Number of indexes on a line
        
        #Curindex = widthPos*lineSize + linePos
        linePos = self.curStep%lineSize #Position along the line
        widthPos = (self.curStep - linePos)/lineSize #Position along the width
        
        lastIndex = self.thickness*lineSize
        
        print linePos, "out of", lineSize        
        print widthPos, "out of", self.thickness       
        
        #Iterate the indexes and keys
        if self.curStep > lastIndex:
            self.curStep = 0
            
            if self.subKeyIndex == len(self.positions[self.keys[self.keyIndex]]) - 1:
                self.subKeyIndex = 0
                self.keyIndex += 1
                
                if self.keyIndex == len(self.keys):
                    self.finished = True
                    
            else:
                self.subKeyIndex += 1
                
        
        pos = nav["init"] + self.separation*(nav["iter"]*linePos + nav["next"]*widthPos)
        self.curStep += 1
        print "trying to place", pos
        
        if self.isOccupied(pos):
            return None
        else:
            self.allPos.append(pos)          
            return pos
        
        
        return array([uniform(0, self.mesh.lx), uniform(0, self.mesh.ly)]) ##DUMMY
            
            

class tableParser():
    
    anyNumber = "\-?\d+\.?\d*[eE]?[+-]?\d*\.?\d*"
    knownSubs = {
                 "lx" : "self.mesh.lx", 
                 "ly" : "self.mesh.ly", 
                 "lz" : "self.mesh.lz"
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
                
                "partition" : ["cyclic"],
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

    def __init__(self, mesh, N):
        
        self.mesh = mesh        
        
        self.N = N        
        
 
    
    
    def unpackMix(self, mixProp):
        
        mixProp = self.getValuesConstraitsEvents(mixProp)
                        
        mixTable = self.initializeTable(mixProp)
        
        
        return mixTable
        
    def initializeTable(self, mixProp):
        
        table = self.unpackStickies(mixProp["fixed"])
        
        
        
        
        return table
    
    
    def unpackStickies(self, stickyProps, appendTo=None):
    
        stickyTable = []
        
        positions = stickyProps["positions"]
        partition = stickyProps["partition"]
        fraction  = stickyProps["fraction"]
        thickness = stickyProps["thickness"]
        separation= stickyProps["separation"]        
        
        nCap = self.N
        
        if fraction != "automatic":
            
            if 1 < fraction or fraction < 0: 
                raise Exception("Fraction must be on [0, 1]")
            elif fraction == 0:
                return []
            
                        
            
            nCap = int(nCap*fraction)
        
        
      
        
        if type(positions) != list:
            if "random" in positions:
                
                try:
                    rngStr = positions.split()[1] 
                except:
                    rngStr = "uniform"          
                    
                try:
                    rng = eval(rngStr)
                except:
                    raise Exception("Unknown rng: " % rngStr)
                
                posIter = ranDist(rng, size=self.mesh.dim) ## Fix input bredde etc.
                    
            else:
                posIter = plcDist(self.mesh, thickness, separation, positions)
        else:
            posIter = plcDist(self.mesh, thickness, separation, *positions)
                
        n = 0
        
        while n < nCap:
            
            newPos = None
            while newPos is None and not posIter.finished:
                newPos = posIter.next()
            
            if posIter.finished:
                break
            
            newVel = normal(size=self.mesh.dim)

            sigma = stickyProps["sigmas"][0]
            eps = stickyProps["epses"][0]
            mass = stickyProps["masses"][0]            
            
            stickyTable.append([False, sigma, eps, mass, newPos, newVel])
      
            n += 1
        print "done"
            
        
        if appendTo:
            return appendTo + stickyTable
        return stickyTable
        
        
        
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