class station:
    def __init__(self,variableNameList,noSpanPts):
        self.names = variableNameList
        self.variables = [[0]*noSpanPts for i in range(len(variableNameList))]
        self.radius = [[0]*noSpanPts]*1
        self.stationType = ""
            
# parse single line of "Axial station" record    
    def loadLine(self,i,string,nameString): # i points at span position
        
        splitLine = nameString.split() # split name String
        valVec = string.split() # split corresponding line of values
        
        for j in range(len(self.names)):
            for k in range(len(splitLine)): 
                if splitLine[k] in self.names[j]:
                    self.variables[j][i] = float(valVec[k]) 
                elif splitLine[k] == "Radius": 
                    self.radius[0][i] = float(valVec[k]) 

# set station type
    def setStationType(self,string): 
        
        if "Duct" in string: 
            self.stationType = "Duct"
        elif "Stator le" in string: 
            self.stationType = "Stator le"
        elif "Stator te" in string: 
            self.stationType = "Stator te"
        elif "Rotor le" in string: 
            self.stationType = "Rotor le"
        elif "Rotor te" in string: 
            self.stationType = "Rotor te"

# get station type
    def getStationType(self):
        return self.stationType


