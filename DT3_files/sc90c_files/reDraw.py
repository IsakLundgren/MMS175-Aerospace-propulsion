import os

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
import station

import sys
import tkinter

path = "./"
file_name = "dt3a.out"


def rePlotResults(variable,full_filename):

  log = False
  stations, noSpanPts, noStations, noBlades = loadStationData(variable,full_filename,log)
  
  if calcFlowFactor: 
    stationVariables = ["flowFactor"]    
    variable = ["flowFactor"]    
    computeFlowFactor(variable,stations,noSpanPts,noStations) # will reduce the list of variables back to a single. Compute flow factors wherever possible.

  if calcStageLoad: 
      stationVariables = ["stageLoad"]    
      variable = ["stageLoad"]    
      computeStageLoad(variable,stations,noSpanPts,noStations) # will reduce the list of variables back to a single. Compute stage loads wherever possible.
  if calcDegreeOfReact: 
    stationVariables = ["degreeOfReact"]    
    variable = ["degreeOfReact"]
    computeDegreeOfReact(variable,stations,noSpanPts,noStations) # will reduce the list of variables back to a single. Compute degree of reaction wherever possible.
  if calcDrosselZiffer: 
    stationVariables = ["DrosselZiffer"]    
    variable = ["DrosselZiffer"]    
    computeDrosselZiffer(variable,stations,noSpanPts,noStations) # will reduce the list of variables back to a single. Compute degree of reaction wherever possible.
  if calcSoverC: 
    stationVariables = ["SoverC"]
    variable = ["SoverC"]
    computeSoverC(variable,stations,noSpanPts,noStations,noBlades) # will reduce the list of variables back to a single. Compute degree of reaction wherever possible.

  showPlots(variable,stations,noSpanPts)

#-------------- sOverc 
def computeSoverC(variable,stations,noSpanPts,noStations,noBlades):

  variable = "s/c"
  for i in range(noStations-1): # load all axial station data
    updateStation0(stations[i],stations[i+1],noSpanPts,noBlades)  
  
  return

def updateStation0(station,stationNext,noSpanPts,noBlades):
  import math
  station.names = "SoverC"
  
  blPtr = 0 
  if "Rotor le" in station.stationType or "Stator le" in station.stationType:
    for i in range(noSpanPts):
      beta1 = station.variables[0][i] ; beta2 = stationNext.variables[0][i]
      r1 = station.variables[1][i] ; r2 = stationNext.variables[1][i]
      O = math.pi*(r1+r2)
      x1 = station.variables[2][i] ; x2 = stationNext.variables[2][i]  
      axialChord = x2 - x1
      ang = math.radians((beta1+beta2)/2.0)
      Chord = axialChord / math.cos(ang)
      station.variables[0][i] = (O/float(noBlades[blPtr]))/Chord
    blPtr = blPtr + 1
      
  return
  
#-------------- flow factor stuff 
def computeFlowFactor(variable,stations,noSpanPts,noStations):

  variable = "flowFactor"
  for i in range(noStations): # load all axial station data
    updateStation1(stations[i],noSpanPts)


def updateStation1(station,noSpanPts):
    station.names = "flowFactor"
    if "Rotor" in station.stationType:
      for i in range(noSpanPts):
        Vm = station.variables[1][i] ; U = station.variables[0][i]
        station.variables[0][i] = station.variables[1][i]/station.variables[0][i]  
        station.variables[1][i] = 0.0  
        if i == getMidPtr(noSpanPts): 
          print("Flow factor at stage mid is = " + str(station.variables[0][i]))
  
#---------------- stage load stuff
def computeStageLoad(variable,stations,noSpanPts,noStations):

  variable = "StageLoad"
  for i in range(noStations-1): # load all axial station data
    updateStation2(stations[i],stations[i+1],noSpanPts)  

def updateStation2(stationLe,stationTe,noSpanPts):
    station.names = "flowFactor"
    if "Rotor le" in stationLe.stationType:
      for i in range(noSpanPts):
        dH = 1000.0*(stationTe.variables[1][i]-stationLe.variables[1][i])
        U = stationLe.variables[0][i]
        stationLe.variables[0][i] = dH/U**2  
        stationLe.variables[1][i] = 0.0  
        if i == getMidPtr(noSpanPts): 
          print("Stage load at stage mid is = " + str(stationLe.variables[0][i]))

#---------------- degree of reaction stuff
def computeDegreeOfReact(variable,stations,noSpanPts,noStations):
  variable = "degreeOfReact"
  for i in range(noStations-1): # load all axial station data
    if "Rotor le" in stations[i].stationType:
      stPos = getStatorLeadingEdgePos(i+2,stations) # search from i+2 and downstream => algorithm will work even if there is a rotor stator intermediate point
      updateStation3(stations[i],stations[i+1],stations[stPos],noSpanPts)  

def getStatorLeadingEdgePos(sPos,stations):
  for i in range(sPos,len(stations)): 
    if "Stator te" in stations[i].stationType: 
      stPos = i
      return stPos

def updateStation3(stationRLe,stationRTe,stationSTe,noSpanPts):
    station.names = "degreeOfReact"
    for i in range(noSpanPts):
      dTRotor = stationRTe.variables[0][i]-stationRLe.variables[0][i]
      dTStage = stationSTe.variables[0][i]-stationRLe.variables[0][i]
      stationRLe.variables[0][i] = dTRotor/dTStage 
      stationRTe.variables[0][i] = 0.0
      stationSTe.variables[0][i] = 0.0
      if i == getMidPtr(noSpanPts):
        print("Degree of reaction at stage mid is = " + str(stationRLe.variables[0][i]))


#-------------- die lokale Drosselziffer berechnen
def computeDrosselZiffer(variable,stations,noSpanPts,noStations):

  variable = "DrosselZiffer"
  for i in range(noStations-1): # load all axial station data
    updateStation4(stations[i],stations[i+1],noSpanPts)  

def updateStation4(stationLe,stationTe,noSpanPts):
    stationLe.names = "DrosselZiffer"; stationTe.names = "DrosselZiffer"
    if "Rotor le" in stationLe.stationType:
      for i in range(noSpanPts):
        U = stationLe.variables[0][i] ; Vm = stationLe.variables[1][i]
        Htot1 = 1000.0*stationLe.variables[2][i] ;  Htot2 = 1000.0*stationTe.variables[2][i]
        psi = (Htot2-Htot1)/U**2 ; fi = Vm/U
        stationLe.variables[0][i] = psi/fi**2  
        stationLe.variables[1][i] = 0.0  ; stationLe.variables[2][i] = 0.0  
        if i == getMidPtr(noSpanPts): 
          print("Flow factor at stage mid is = " + str(stationLe.variables[0][i]))
#----

def getMidPtr(noSpanPts):    
    return int((noSpanPts-1)/2)

def loadStationData(variable,ffn,log): # ffn = full file name
  
  try:
    file = open(ffn,mode='r')
  except IOError:
    print("Failed to open file"+ffn)
    
# get number of stations and streamlines. Reposition to first line 
  noStations, noSpanPts = getDimensions(file,log)
   
  stationList = [] # init a list
  for i in range(noStations): # load all axial station data
    stationList.append(station.station(variable,noSpanPts)) # create and initialize a new station in the list of objects
    loadRecord(i,file,stationList[i],noSpanPts) # load record i from file to recordList by parsing all records using parseLine     

  loadLosses(file,stationList,noSpanPts) # load losses from file      
  
  file.seek(0)  
  noBlades = getNoBlades(file,stationList,noSpanPts) # load losses from file        

  file.close()
  
  return (stationList, noSpanPts, noStations, noBlades)

def countNumberOfPlots(variable,stations,plotSelection):
  
  countNumberOfPlots = 0 
  for station in stations: 
    if plotShouldBeIncluded(variable,station.stationType,plotSelection): 
      countNumberOfPlots = countNumberOfPlots + 1
    
  return countNumberOfPlots

def plotShouldBeIncluded(variable,sType,plotSelection):

  plotShouldBeIncluded = False
  
  if plotSelection == "All":
    plotShouldBeIncluded = True
  elif plotSelection == "Blades": 
    if sType =="Stator le" or sType =="Stator te" or \
       sType=="Rotor le" or sType =="Rotor te":
      plotShouldBeIncluded = True
  elif plotSelection == "TrailingEdges": 
    if sType=="Stator te" or sType =="Rotor te":
      plotShouldBeIncluded = True      
  elif plotSelection == "LeadingEdges": 
    if sType=="Stator le" or sType =="Rotor le":
      plotShouldBeIncluded = True      
  elif plotSelection == "Rotors":
    if sType=="Rotor le" or sType =="Rotor te":
      plotShouldBeIncluded = True
  elif plotSelection == "Stators": 
    if sType =="Stator le" or sType =="Stator te":
      plotShouldBeIncluded = True
  elif plotSelection == "RotorTrailingEdges": 
    if sType =="Rotor te":
      plotShouldBeIncluded = True      
  elif plotSelection == "RotorLeadingEdges": 
    if sType =="Rotor le":
      plotShouldBeIncluded = True      
      
  return plotShouldBeIncluded

  
def getSubPlotDistribution(noPlots):
  import math
  
  noRows = int(math.ceil(math.sqrt(noPlots)))
  noColumns = noRows
        
  return (noRows, noColumns)

def getBLNam(noStators,noRotors,stationType): 

# Copy by default
  getBLNam = stationType
  
  if "Stator" in stationType:
    if noRotors == 0: 
      getBLNam = "IGV "
    else:
      BladePos = stationType.split()[1]
      getBLNam = "Stator "+str(noStators)+" "+BladePos
      
  if "Rotor" in stationType:
    BladePos = stationType.split()[1]
    getBLNam = "Rotor "+str(noRotors)+" "+BladePos
      
  return getBLNam

def showPlots(variable,stations,noSpanPts):

  stationVariables, plotSelection = getRunOptions("runOptions.txt")
  noPlots = countNumberOfPlots(variable,stations,plotSelection)
  noRows, noColumns = getSubPlotDistribution(noPlots)

  counter = 0 
  noStatorRows = 0 
  noRotorRows = 0 
  plotPtr = 0 
  
  for i in range(len(stations)): 
# keep track of number of stator rows / rotor rows (used to differentiate between IGV and stator
    if "Stator le" in stations[i].stationType:
      if noRotorRows > 0: # this cant be IGV
        noStatorRows = noStatorRows + 1 
    if "Rotor le" in stations[i].stationType: 
      noRotorRows = noRotorRows + 1 
#
    if plotShouldBeIncluded(variable,stations[i].stationType,plotSelection):
      plt.subplot(noRows, noColumns, int(counter)+1)
      plt.plot(stations[i].variables[plotPtr][0:noSpanPts], \
               stations[i].radius[0:noSpanPts][0],'bo-')
      counter = counter + 1
      plt.tick_params(labelsize=8)
      plt.subplots_adjust(wspace=0.35, hspace=0.35)
      plt.grid(True)
      titleStr1 = getBLNam(noStatorRows,noRotorRows,stations[i].stationType)  
      plt.title(titleStr1)
  
##  mng = plt.get_current_fig_manager()
##  mng.window.state('zoomed') #works fine on Windows!

  
  plt.suptitle(variable[0], fontsize=36)

##  plt.show()
  print("Completed plotting")
  

def getDimensions(file,log):
  noStations = 0 
  noSpanPts = 0 
  dataSet = False
  
  if log: 
    print("attempting to retrieve dimensions. opening file",file)
  
  while not dataSet:      
    try: # read until error occurs
      string = file.readline().rstrip() # to avoid breaking on an empty line
    except IOError:
      break
# stations
    if "Ax dist hub" in string: # parse out number of stations
        if log: 
          print("found ax dist hub location")
        next(file) # skip empty line
        eos = False # end of stations
        while not eos:
          string = file.readline().rstrip()
          if string =="":
            eos = True
          else:
            noStations = int(string.split()[0])
# span positions
    if "Rad Stn" in string: # parse out number of stations
      if log: 
        print("rad stn location")
      next(file) # skip empty line
      next(file) # skip empty line
      eos = False # end of stations
      while not eos:
        string = file.readline().rstrip()
        if string =="":
          eos = True
        elif "mass avge" in string:
          break
        else:
          noSpanPts = int(string.split()[0])   
    
    if noStations > 0 and noSpanPts > 0: 
      dataSet = True

  if log:
    print("file parsing complete")
  file.seek(0)
  
  return (noStations, noSpanPts)

def isALossVariable(varName):
  
  isALossVariable = False
  if varName == ["total"]:
    isALossVariable = True
  elif varName == ["profile"]:
    isALossVariable = True
  elif varName == ["shock"]:
    isALossVariable = True
  elif varName == ["secondary"]:  
    isALossVariable = True
  elif varName == ["tip-clearance"]:  
    isALossVariable = True
  
  return isALossVariable
  
def loadRecord(i,file,station,nSP): 
  
  while 1:      
    try: # read until error occurs
      string = file.readline().rstrip() # to avoid breaking on an empty line
    except IOError:
      break
# stations
    if string[0:13]=="Axial station":
      parseStationData(i,file,station,nSP)
      break
    elif "Breakdown of loss coefficients" in string: 
      break # last time in loop
  
  return 

def loadLosses(file,stations,nSP): 
  
  while 1:      
    try: # read until error occurs
      string = file.readline().rstrip() # to avoid breaking on an empty line
    except IOError:
      break
 
    if string[0:6]==" Plane":
      parseProfileData(file,stations[getIndex(string)],nSP)
    elif "SUMMARY OF STAGE PERFORMANCE" in string: 
      break
    
  return 


def getNoBlades(file,stations,nSP): 
  
  noBlades = []
  while 1:      
    try: # read until error occurs
      string = file.readline().rstrip() # to avoid breaking on an empty line
    except IOError:
      break
 
    if "NO.OF BLADES" in string:
      noBlades.append(string.split()[3])
    elif "Axial station" in string:  
      break
    
  return noBlades

def getIndex(string): 
  getIndex = int(string.split()[1])-1
  
  return getIndex

def parseLossBlock(file,station,nameString,nSP):
  next(file) # skip lines to first set of data

  for i in range(nSP):
    string = file.readline().rstrip() 
    station.loadLine(i,string,nameString)


def parseProfileData(file,station,nSP):

  next(file)
  string = file.readline().rstrip() # to avoid breaking on an empty line
  
  parseLossBlock(file,station,mergeWords(string),nSP)
        
def parseStationData(i,file,station,nSP):

  string = file.readline().rstrip() # to avoid breaking on an empty line
  station.setStationType(string)
  noParsedBlocks = 0
  totNoBlocks = getNoBlStType(station.getStationType()) 
  
  print(station.stationType+" is at plane "+str(i+1))
  
  while 1:
    string = file.readline().rstrip() # to avoid breaking on an empty line
    if "Rad Stn" in string or "Mcrit" in string or "rVswirl" in string or \
          "Maxial" in string or string[0:5]==" Mabs": # parse out number of stations (first block of data)
      parseBlock(file,station,mergeWords(string),nSP)
      noParsedBlocks = noParsedBlocks + 1 
      if noParsedBlocks == totNoBlocks:
        break

def getNoBlStType(stype):
  
  if stype == "Duct":
    getNoBlStType = 3
  elif stype == "Stator le":
    getNoBlStType = 3    
  elif stype == "Stator te":
    getNoBlStType = 4
  elif stype == "Rotor le":
    getNoBlStType = 4    
  elif stype == "Rotor te":
    getNoBlStType = 5
  
  return getNoBlStType

      
def mergeWords(string):  # replace selected spaces with hyphen => parser does not confuse the number of words in the line
  mergeWords = string

# fix space in "Rad Stn"
  if "Rad Stn" in string: 
    pos = string.find("Rad Stn")
    mergeWords = string[0:pos+3] + "-" + string[pos+4:len(string)]
  
  if "tip clearance" in string: 
    pos = string.find("tip clearance")
    mergeWords = string[0:pos+3] + "-" + string[pos+4:len(string)]
  
  if "Ax dist" in string: 
    pos = mergeWords.find("Ax dist")  
    mergeWords = mergeWords[0:pos+2] + "-" + mergeWords[pos+3:len(string)]
  
  if "de Haller" in string: 
    pos = mergeWords.find("de Haller")  
    mergeWords = mergeWords[0:pos+2] + "-" + mergeWords[pos+3:len(string)]

  if "ro x Vm" in string: 
    pos = mergeWords.find("ro x Vm")  
    mergeWords = mergeWords[0:pos+2] + "-" + mergeWords[pos+3:len(string)]
    mergeWords = mergeWords[0:pos+4] + "-" + mergeWords[pos+5:len(string)]

  if "ax force" in string: 
    pos = mergeWords.find("ax force")  
    mergeWords = mergeWords[0:pos+2] + "-" + mergeWords[pos+3:len(string)]

  if "tang force" in string: 
    pos = mergeWords.find("tang force")  
    mergeWords = mergeWords[0:pos+4] + "-" + mergeWords[pos+5:len(string)]

  if "blade angle" in string: 
    pos = mergeWords.find("blade angle")  
    mergeWords = mergeWords[0:pos+5] + "-" + mergeWords[pos+6:len(string)]

  if "del H" in string: 
    pos = string.find("del H")
    mergeWords = mergeWords[0:pos+3] + "-" + mergeWords[pos+4:len(string)]
    
  return mergeWords
  
def parseBlock(file,station,nameString,nSP):
  next(file) ; next(file) # skip lines to first set of data

  for i in range(nSP):
    string = file.readline().rstrip() 
    station.loadLine(i,string,nameString)

def getStatus(path):
  try:
    file = open(path+"status.txt",mode='r')
  except OSError: 
    print("failed to open file")
    exit(0)
    
  getStatus = int(file.readline().rstrip())
  file.close()
  
  return getStatus # to avoid breaking on an empty line

def setStatus(istat,path):
  try:
    file = open(path+"status.txt",mode='w')
  except OSError: 
    print("failed to open file")
    exit(0)    
    
  file.write(str(istat)+"\n") 
  file.close()

#def sel():
#  global stationVariables
#  stationVariables = [MODES[int(v.get())]]
#  rePlotResults(stationVariables,path+file_name)
#
#def sel2():
#  global plotSelection
#  plotSelection = TYPES[int(v.get())]
#  rePlotResults(stationVariables,path+file_name)
#
#
#
#
#loadOptions()
def getRunOptions(fn):
  global calcFlowFactor
  global calcStageLoad
  global calcDegreeOfReact
  global calcDrosselZiffer
  global calcSoverC

  calcFlowFactor = False
  calcStageLoad = False
  calcDegreeOfReact = False
  calcDrosselZiffer = False
  calcSoverC = False

  try:
    file = open(path+fn,mode='r')
  except OSError: 
    print("failed to open file"+fn)
    exit(0)    
    
  EOF = False
  while not EOF: 
    try: # read until error occurs
      string = file.readline().rstrip() # to avoid breaking on an empty line
    except EOFError:
      break
    
    if "stationVariableSelection:" in string:
      pos = string.find(":")
      if pos > 0: 
        stationVariables = [string[pos+1:len(string)].lstrip()]
# special handling (search for this tag to extend functionality)
      if stationVariables[0] == "flowFactor": 
        stationVariables[0] = "U"
        stationVariables.append("Vm")
        calcFlowFactor = True
      if stationVariables[0] == "stageLoading":
        stationVariables[0] = "U"
        stationVariables.append("Htot")
        calcStageLoad = True
      if stationVariables[0] == "degreeOfReaction":
        stationVariables[0] = "Tstat"
        calcDegreeOfReact = True
      if stationVariables[0] == "DrosselZiffer":
        stationVariables[0] = "U"
        stationVariables.append("Vm")
        stationVariables.append("Htot")
        calcDrosselZiffer = True
      if stationVariables[0] == "SoverC":
        stationVariables[0] = "blade-angle"
        stationVariables.append("Radius")
        stationVariables.append("Ax-dist")
        calcSoverC = True
    
    if "plotSelection:" in string:
      pos = string.find(":")
      if pos > 0: 
        plotSelection = string[pos+1:len(string)].lstrip()
    
    if string == "": 
      EOF = True
      
  if len(stationVariables[0].split()) == 2: 
    stationVariables[0] = \
      stationVariables[0].split()[0]+"-"+stationVariables[0].split()[1]

  return stationVariables, plotSelection

#
# Main program
#

##stationVariables, plotSelection = getRunOptions("runOptions.txt")
##rePlotResults(stationVariables,path+file_name)
