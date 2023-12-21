#
# 2014-12-17 Program that generates the indata for the SC90C (MATLAB).
#            Programmer: Egill Thorbergsson. 
# 
# 2015-02-22 Program transferred to and tailored for python environment. 
#            Programmer: Tomas Gronstedt
#
import sc90c_simulation_input
import sys, os
from subprocess import call
import numpy
import station
import matplotlib
matplotlib.use('TkAgg')

from shutil import copyfile
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import tkinter
from tkinter import *
from tkinter import ttk

name = "dt3a" # base name for SC90C files
path = "./" # path
inpFileName = "smoothInput.txt"
inpFileNameBeforeOpt = "smoothInputBeforeOpt.txt"
inpFileNameOriginal = "smoothInputOriginal.txt"
copyfile(path + inpFileName, path + inpFileNameBeforeOpt)
copyfile(path + inpFileName, path + inpFileNameOriginal)

#### Variables for trackGUI.py###
# path = "sc90c_files/"
file_name = "runOptions.txt"
stationVariables = "Vtot" # init plotting variable choice
plotSelection = "Blades"  # init plotting station choice
lossVariables = "Total" # init plotting loss choice
displayMassAvgVal = False # display mass averaged variable 
TYPES = [("All"), ("Blades"),("LeadingEdges"),("TrailingEdges"),("Rotors"),("Stators"),("RotorLeadingEdges"),("RotorTrailingEdges")]
VARIABLES = [("Vtot"), ("Vm"),("Vrad"),("Vswirl"),("Mabs"),("Vax"),("Maxial"),("Ang.sw"),
         ("Ang.rad"),("Ptot"),("Pstat"),("Ttot"),("Tstat"),("rVswirl"),("Devn."),
         ("Inc"),("U"),("Prel"),("Ang.rel"),("Mrel"),("Vrel"),("de Haller"),("SoverC"),
         ("diffusion"),("flowFactor"),("stageLoading"),("degreeOfReaction"),("Inc")]
LOSSES = [("total"),("profile"),("shock"),("secondary"),("tip-clearance")]

#### End ###


prat = 0.0 ; NH = 0.0 ; mdot = 0.0 ; T0 = 0.0 ; P0 = 0.0
reDistFact  = 0.0 ; n_stages = 0 ; n_span = 0
fl_igv = 0 ; n_bezier = 0 ; n_duct_inlet = 0  ; n_duct_outlet = 0 
n_blade_intermediate_stations = 0 ; hubStretch = [] 
tipStretch = [] ; stretch = []
blades = []
xbezHub = [] ; rbezHub = [] ; xbezShroud = [] ; rbezShroud = []
spacing = [] ; hubSpacing = [] ; tipSpacing = []
workFact_R1 = []
alpha_TE_S1 = []
alpha_TE_IGV = []
xTip = [] ; xHub = []
inpData = [("prat:",prat),
           ("NH:",NH),
           ("mdot:",mdot),
           ("T0:",T0),
           ("P0:",P0),
           ("reDistFact:",reDistFact),
           ("n_stages:",n_stages),
           ("n_span:", n_span),
           ("fl_igv:", fl_igv),
           ("n_bezier:", n_bezier),           
           ("n_duct_inlet:", n_duct_inlet),
           ("n_duct_outlet:", n_duct_outlet),
           ("n_blade_intermediate_stations:", n_blade_intermediate_stations),
           ("hubStretch:", hubStretch),
           ("tipStretch:", tipStretch),
           ("stretch:", stretch),
           ("blades:", blades),
           ("xbezHub:",xbezHub),
           ("rbezHub:",rbezHub),
           ("xbezShroud:",xbezShroud),
           ("rbezShroud:",rbezShroud),
           ("spacing:",spacing),
           ("hubSpacing:",hubSpacing),
           ("tipSpacing:",tipSpacing),
           ("alpha_TE_IGV:",alpha_TE_IGV),
           ("alpha_TE_S1:",alpha_TE_S1),
           ("workFact_R1:",workFact_R1),
           ("xTip:",xTip),
           ("xHub:",xHub)
           ]

xmin = 0.0 ; xmax = 0.0
sVecHubwWB = [] ; eVecHubwWB = [] 
sVecTipwWB = [] ; eVecTipwWB = [] 


streamLines = []

def setGasPathPositions(compressor,deltaTot,xHub,xTip,n_duct_inlet,n_duct_outlet):

    xLength = xHub[-1]-xHub[0] # length of compressor

# set inlet station data
    dxInlet = 0.12*xLength # use 12% of total compressor length to distribute inlet stations
    xInlet = [0.0]*n_duct_inlet
    for i in range(n_duct_inlet):
        iInv = n_duct_inlet-i
        xVal = xHub[0] - dxInlet*float(iInv)/float(n_duct_inlet)
        xInlet[i] = xVal    
    compressor.setDuctInlet(n_duct_inlet,xInlet) # inlet duct    
    xmin = xInlet[0]

# set station data in bladed region
    compressor.setDuctBlades(xHub,xTip) 
    
# set outlet station data
    dxOutlet = 0.12*xLength # use 12% of total compressor length to distribute outlet stations
    xOutletHub = [0.0]*n_duct_outlet
    xOutletTip = [0.0]*n_duct_outlet
    for i in range(n_duct_outlet):
        xVal = xHub[-1] + dxInlet*float(i+1)/float(n_duct_outlet)
        xOutletHub[i] = xVal    
        xVal = xTip[-1] + dxInlet*float(i+1)/float(n_duct_outlet)
        xOutletTip[i] = xVal    
    compressor.setDuctOutlet(n_duct_outlet,xOutletHub,xOutletTip) # outlet duct    
    xmax = max(xOutletTip[n_duct_outlet-1],xOutletHub[n_duct_outlet-1])
    
    return (xmin,xmax)


def initCompressor(path,name):       
    global xTip, xHub, xmin, xmax, compressor
    
    
# construct compressor object.
    compressor = sc90c_simulation_input.sc90cInput(path,name,getInpDataVal("n_stages:"),getInpDataVal("fl_igv:"),getInpDataVal("n_span:"), \
                   getInpDataVal("n_duct_inlet:"),getInpDataVal("n_duct_outlet:"),getInpDataVal("n_blade_intermediate_stations:"),getInpDataVal("n_bezier:")) # creates ISIG through constructor method

    xTip = getInpDataVec("xTip:") ; stretch = getInpDataVec("stretch:") ; spacing = getInpDataVec("spacing:")
    
    xTip, deltaTot = compressor.stretchChord(xTip,stretch,spacing) # stretch original chord and spacing at mid
    nrows = 2*(getInpDataVal("fl_igv:") + 2*getInpDataVal("n_stages:")) 
    xHub[0:nrows]= xTip[0:nrows]
    
    xHub, xTip = applyHubAndTipStretch(getInpDataVec("hubStretch:"),getInpDataVec("tipStretch:"),xHub,xTip)    
    xHub, xTip = applyHubAndTipSpacing(getInpDataVec("hubSpacing:"),getInpDataVec("tipSpacing:"),xHub,xTip) # allows sweep
    
    xmin, xmax = setGasPathPositions(compressor,deltaTot,xHub,xTip,getInpDataVal("n_duct_inlet:"),getInpDataVal("n_duct_outlet:"))
        
# init compressor     
    compressor.Rgas = 0
    
    compressor.setOptVars(getInpDataVal("prat:"),getInpDataVal("reDistFact:")) # set optimization data
    compressor.setPerformance(getInpDataVal("NH:"),getInpDataVal("mdot:"),getInpDataVal("T0:"),getInpDataVal("P0:"))
    compressor.setPressRatios(getInpDataVal("prat:"),getInpDataVal("reDistFact:"),getInpDataVec("workFact_R1:")) # PR_R1    
    compressor.setSwirl(getInpDataVec("alpha_TE_IGV:"),getInpDataVec("alpha_TE_S1:"))
    
    compressor.setNoBlades(getInpDataVec("blades:")) # IGV, R1, S1
    # the function below should update the xbez to match the desired corner points
    compressor.set_Bezier_x_Points(xHub,xTip) # Bezier curve control points (axial coordinate)
    compressor.set_Bezier_r_CurvePoints(getInpDataVec("rbezHub:"),getInpDataVec("rbezShroud:")) # Bezier curve control points (radial coordinate)
            
    return compressor

def getXBezier(xmin,xIGVTe,xR1Te,xmax): # control using same x-coordinates for hub and shroud
    # first point at entry
    # second point at IGV exit (since hub is intended to have an inflexion point there)
    # remaining points distributed linearly

    nBez = getInpDataVal("n_bezier:")
    
    if nBez == getNoXInFullGrid(): # inlet + outlet + one behind each blade row
        # interpolation to full grid already happened
        return getFullXBezier(xmin,xmax,xHub,xTip)
    else:
        xBezier = [0.0]*nBez
        xBezier[0] = xmin
        xBezier[1] = xIGVTe
        xBezier[2] = xIGVTe
        
        xpatch = numpy.linspace(xIGVTe,xmax,nBez-2).tolist()
        xBezier[3:nBez] = xpatch[1:nBez-2]
        
        return (xBezier,xBezier)

# function determines how many x-points there should be when there is max. control on the Bezier curve 
def getNoXInFullGrid(): 
    nox = 2 + (getInpDataVal("fl_igv:") + 2*getInpDataVal("n_stages:") - 1) + 2
    
    return nox


def applyHubAndTipStretch(hStretch,tStretch,xHub,xTip): 
    
    counter = 0
    for i in range(len(xTip)): 
        if i % 2 == 0: # if i is even => leading edge. (0,2,4,6,8,... are leading edges)
            xHubLe = xHub[i] ; xHubTe = xHub[i+1]
            xmid = (xHubLe + xHubTe)/2.0
            xHub[i]   =   xHub[i] + (hStretch[counter]-1.0)*(xHub[i]   - xmid)
            xHub[i+1] = xHub[i+1] + (hStretch[counter]-1.0)*(xHub[i+1] - xmid)

            xTipLe = xTip[i] ; xTipTe = xTip[i+1]
            xmid = (xTipLe + xTipTe)/2.0
            xTip[i]   =   xTip[i] + (tStretch[counter]-1.0)*(xTip[i]   - xmid)
            xTip[i+1] = xTip[i+1] + (tStretch[counter]-1.0)*(xTip[i+1] - xmid)

            counter = counter + 1
    
    return (xHub, xTip)


def applyHubAndTipSpacing(hSpacing,tSpacing,xHub,xTip): 

    # if all elements in hSpacing and tSpacing are one, xHub/xTip are left unaffected
    bladeCounter = 0 ; noX = len(xTip) ; deltaTot = 0.0
    for i in range(2,noX): # first two are first blade => no upstream spacing.
        if i % 2 == 0: # if i is even => leading edge. (0,2,4,6,8,... are leading edges)
            delta = (xHub[i] - xHub[i-1])*(hSpacing[bladeCounter]-1.0)
            for j in range(i,noX): # add delta to all downstream x-coordinates at Hub
                xHub[j] = xHub[j] + delta
            deltaTot = deltaTot + delta
            bladeCounter = bladeCounter + 1
            
    # if all elements in hSpacing and tSpacing are one, xHub/xTip are left unaffected
    bladeCounter = 0 ; noX = len(xTip) ; deltaTot = 0.0
    for i in range(2,noX): # first two are first blade => no upstream spacing
        if i % 2 == 0: # if i is even => leading edge. (0,2,4,6,8,... are leading edges)
            delta = (xTip[i] - xTip[i-1])*(tSpacing[bladeCounter]-1.0)
            for j in range(i,noX): # add delta to all downstream x-coordinates at Hub
                xTip[j] = xTip[j] + delta
            deltaTot = deltaTot + delta
            bladeCounter = bladeCounter + 1
    
    
    return (xHub, xTip)

def executeSC90C(): 
    import win32com.client 
    import time

    shell = win32com.client.Dispatch("WScript.Shell")
    shell.AppActivate("Command Prompt")
    
    shell.SendKeys("runSC90C.cmd")
    shell.SendKeys("{ENTER}")
    time.sleep(0.5) # delays for 0.5 seconds. This avoids having keys sent to debug environment & reading results before the SC90C code has finished
    
    shell.SendKeys("y")    
    shell.SendKeys("{ENTER}")
    shell.SendKeys("{ENTER}")
    
def loadData(fn): 
    try:
        file = open(path+fn,mode='r')
    except IOError:
        print("Failed to open file"+path+fn)
    while 1: 
        string = file.readline().lstrip() # to avoid breaking on an empty line
        if string == "":
            break
        elif string[0] == "#": 
            continue
        else:
            loadInputVariable(string,inpData)
    
    file.close()


def loadInputVariable(string,inpData):
    
    counter = 0
    for name,variable in inpData:
        if name in string:
            splitData = string.split()
            element = list(inpData[counter]) # cast tuple to a list
            element[1] = splitData[1:len(splitData)]
            inpData[counter] = tuple(element)
            break
        counter = counter + 1
            
    return 

#
# Main program: 
#
#         * Only data that will be part of optimization is visible in main program.
#
def main(includePlotting): 
    global streamLines



# returns a compressor object    
    compressor = initCompressor(path,name) # all data initialized to values
    
    compressor.calculation()
    
    compressor.createSC90C_Input_Files()
##    global ars
##    global arr
    ars, arr = compressor.logAspectRatios()
    # compressor.writeCFDInput()
        
    executeSC90C() 


    retrieveValues()            
 
    if includePlotting:
        f = plt.figure(1)
        compressor.createPlots(streamLines) # create gas path plot (based on input)
        canvas = FigureCanvasTkAgg(f, master=CompFigure)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky='EWNS',padx=5, pady=30)

    ###---------- Summary of the computation -------------###
    ComputationSummary=LabelFrame(masterPost, text="Summary of the computation")
    ComputationSummary.grid(column=1, row=0, sticky='NEWS')

    file = open(path+name+".out",mode='r')

    if not execSuccess(file): 
        ErrorTitle = Label(ComputationSummary, text="The case has not converged, the simulation has failed. Try another geometry." )
        ErrorTitle.grid(row=0, column=0, sticky='NWS', padx=15)
##        Error.grid(row=0, column=1, sticky='NWS', padx=15)
        
    elif "MASS FLOW WAS REDUCED TO" in open(path+name+".out").read():
        ErrorTitle = Label(ComputationSummary, text="The compressor stage is choked. This will not work in CFX. Try another geometry." )
        ErrorTitle.grid(row=0, column=0, sticky='NWS', padx=15)
##        Error.grid(row=0, column=1, sticky='NWS', padx=15)
        
    else:

        StatAspRatTitle = Label(ComputationSummary, text="The aspect ratios are:  " )
        StatAspRat=Label(ComputationSummary, text='IGV:'+ str(round(ars[0],2))+'  R1:'+ str(round(arr[0],2)) +'  S1:'+ str(round(ars[1],2)) )
        StatAspRatTitle.grid(row=0, column=0, sticky='NWS', padx=15)
        StatAspRat.grid(row=0, column=1, sticky='NWS', padx=15)


        StageEfficiencyTitle =  Label(ComputationSummary, text="Stage efficiencies are:")
        StageEfficiencyTitle.grid(row=1, column=0, sticky='NWS', padx=15)
        StageEfficiency =  Label(ComputationSummary, text=str(stageEff[0]))
        StageEfficiency.grid(row=1, column=1, sticky='NWS', padx=15)
        
        EfficiencyPolTitle = Label(ComputationSummary, text="The computed stage efficiency in SC90C is:")
        EfficiencyPolTitle.grid(row=2, column=0, sticky='NWS', padx=15)
        EfficiencyPol = Label(ComputationSummary, text=str(effPol))
        EfficiencyPol.grid(row=2, column=1, sticky='NWS', padx=15)



def updateModel():

    readInputData(inpFileName)

    #UpdateModelJoint
    nBlRows = getInpDataVal("fl_igv:") + 2*getInpDataVal("n_stages:")         

    for i in range(nBlRows): 
        setInpData("stretch:",eVecwWC[i].get(),i)
        # Update position of xcoordinates (same as new corner points)
        # x1 - inlet; x2 - igv_te; x3 - r1_te; x4 - outlet
        #setInpData("xbezHub:", [updated vector], i)
        #setInpData("xbezShroud:", [updated vector], i)
  
    for i in range(nBlRows): 
        setInpData("blades:",eVecwWB[i].get(),i)        
  
    for i in range(nBlRows): 
        setInpData("hubStretch:",eVecHubwWHS[i].get(),i)        
        setInpData("tipStretch:",eVecTipwWHS[i].get(),i)        

    for i in range(nBlRows): 
        setInpData("spacing:",eVecwWS[i].get(),i)        

    #UpdateBez

    for i in range(getInpDataVal("n_bezier:")): 
        setInpData("rbezHub:",eVecHubwWB[i].get(),i)        
        setInpData("rbezShroud:",eVecTipwWB[i].get(),i)        

    #UpdateModelStage
    setInpData("alpha_TE_IGV:", eVecHubwWS[0].get(), 0);
    setInpData("alpha_TE_IGV:", eVecMidwWS[0].get(), 1);
    setInpData("alpha_TE_IGV:", eVecTipwWS[0].get(), 2)
    setInpData("alpha_TE_S1:", eVecHubwWS[1].get(),0)
    setInpData("alpha_TE_S1:", eVecMidwWS[1].get(),1)
    setInpData("alpha_TE_S1:", eVecTipwWS[1].get(),2)

    #UpdateModelWork
    
    setInpData("workFact_R1:", eVecHubwWW[0].get(),0)
    setInpData("workFact_R1:",eVecMidwWW[0].get(),1)
    setInpData("workFact_R1:", eVecTipwWW[0].get(),2)

    main(True)
    
def updateAirFile(): 
# read aor-file
    fn = path+name+".aor" ; aorFile = []
    try:
        file = open(fn,mode='r')
    except IOError:
        print("Failed to open file"+fn)
    
    while 1: 
        string = file.readline().rstrip() 
        if string == "":
            break   
        else: 
            aorFile.append(string + " \n") # to avoid breaking on an empty line            
    
    file.close()
    
    return

def updateGeoFile():
# create copy of file.....
    fn = path+name+".geo" ; geoFile = []
    try:
        file = open(fn,mode='r')
    except IOError:
        print("Failed to open file"+fn)
    
    while 1: 
        string = file.readline().rstrip() 
        if string == "":
            break   
        else: 
            geoFile.append(string + " \n") # to avoid breaking on an empty line            
    file.close()

# replace first line and re-write
    geoFile[0] = " 5   \n" 
    
    fn = path+name+".geo" 
    try:
        file = open(fn,mode='w')
    except IOError:
        print("Failed to open file"+fn)

    for i in range(len(geoFile)): 
        file.write(geoFile[i])

        
    file.close()
    
    return
    
def retrieveValues():
    
    log = False
    try:
        file = open(path+name+".out",mode='r')
    except IOError:
        print("Failed to open file"+path+name+".out")
    
    if not execSuccess(file): 
        print("Calculation failed")
        return

# get number of stations and streamlines. Reposition to first line 
    noStations, noSpanPts = getDimensions(file,log)

    variable = ["Mrel"]
    
    stationList = [] # init a list
    for i in range(noStations): # load all axial station data
        stationList.append(station.station(variable,noSpanPts)) # create and initialize a new station in the list of objects
        loadRecord(i,file,stationList[i],noSpanPts) # load record i from file to recordList by parsing all records using parseLine     
   
    loadLosses(file,stationList,noSpanPts) # load losses from file  
    
    loadOverallPerformance(file,stationList,noSpanPts)
    
      
# parse from top of file again.
    file.seek(0)
    assembleStreamLines(file,noStations,noSpanPts) # load losses from file    

    file.close()

    return 


def computeDrosselZiffer(file,noStations,noSpanPts): # load losses from file

    localStationList = [] # init a list
    for i in range(noStations): # load all axial station data
        localStationList.append(station.station( [("Vm"), ("U"),("Htot")],noSpanPts))
        loadRecord(i,file,localStationList[i],noSpanPts) # load record i from file to recordList by parsing all records using parseLine

    drosselZiffer = getDrosselZiffer(localStationList,noStations,noSpanPts)
    
    return drosselZiffer

def assembleStreamLines(file,noStations,noSpanPts): # load losses from file
    global streamLines

    streamLines = [] # init a list
    for i in range(noStations): # load all axial station data
        streamLines.append(station.station( [("Ax-dist"), ("Radius")],noSpanPts))
        loadRecord(i,file,streamLines[i],noSpanPts) # load record i from file to recordList by parsing all records using parseLine

    return

    
def getDrosselZiffer(localStationList,noStations,noSpanPts): 

    for i in range(noStations): # load all axial station data
        if "Rotor le" in localStationList[i].stationType:
            firstRotorPtr = i
            break

    for i in range(firstRotorPtr+1,noStations): # load all axial station data
        if "Rotor le" in localStationList[i].stationType:
            lastRotorPtr = i
    
    midPtr = getMidPtr(noSpanPts)
            
# get variables at inflow (Grieb 5.2.2.2)
    fiE = localStationList[firstRotorPtr].variables[0][midPtr]/localStationList[firstRotorPtr].variables[1][midPtr]
    fiA = localStationList[lastRotorPtr+1].variables[0][midPtr]/ localStationList[lastRotorPtr+1].variables[1][midPtr]
    fi = (fiE+fiA)/2.0
    
    h02 = localStationList[lastRotorPtr+1].variables[2][midPtr] ; h01 = localStationList[firstRotorPtr].variables[2][midPtr]
    Heff = h02 - h01
    
    U2Sum = 0.0
    for i in range(noStations): # load all axial station data
        if "Rotor le" in localStationList[i].stationType:
            U = localStationList[i].variables[1][midPtr]
            U2Sum = U2Sum + U**2
    
    psi = 1000.0*(2.0*Heff)/U2Sum
    
    
    return (psi/fi**2)
   

def getMidPtr(noSpanPts):
    
    return int((noSpanPts-1)/2)

def loadOverallPerformance(file,stations,nSP):
    global effPol
    global stageEff
    next(file) ; next(file)    

    stageEff = [0.0]*0 ; stageCounter = 0 
    while 1: 
        string = file.readline().rstrip() # to avoid breaking on an empty line
        splitData = string.split()
        if len(splitData) > 3: 
            stageEff.append(float(splitData[3]))
            stageCounter = stageCounter + 1
        else: 
            break

    next(file) ; next(file)    
    string = file.readline().rstrip() # to avoid breaking on an empty line
    
    effPol = float(string.split()[5])

    
# print overall efficiency to standard output. 
    sys.stdout.write('Stage efficiencies are: ')
    for i in range(stageCounter):         
        sys.stdout.write("{:8.4f}".format(stageEff[i]))
        
    print("\n")
    print("Polytropic compressor efficiency",effPol)


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

def getIndex(string): 
    getIndex = int(string.split()[1])-1
  
    return getIndex

def parseProfileData(file,station,nSP):
    
    next(file)
    string = file.readline().rstrip() # to avoid breaking on an empty line
  
    parseLossBlock(file,station,mergeWords(string),nSP)

def parseLossBlock(file,station,nameString,nSP):
    next(file) # skip lines to first set of data
    
    for i in range(nSP):
        string = file.readline().rstrip() 
        station.loadLine(i,string,nameString)


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

def parseStationData(i,file,station,nSP):
    
    string = file.readline().rstrip() # to avoid breaking on an empty line
    station.setStationType(string)
    noParsedBlocks = 0
    totNoBlocks = getNoBlStType(station.getStationType()) 
    
    while 1:
        string = file.readline().rstrip() # to avoid breaking on an empty line
        if "Rad Stn" in string or "Mcrit" in string or "rVswirl" in string or \
          "Maxial" in string or string[0:5]==" Mabs": # parse out number of stations (first block of data)
            parseBlock(file,station,mergeWords(string),nSP)
            noParsedBlocks = noParsedBlocks + 1 
            if noParsedBlocks == totNoBlocks:
                break

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

def execSuccess(file):
    
    execSuccess = True
    while 1:      
        try: # read until error occurs
            string = file.readline().rstrip() # to avoid breaking on an empty line
        except IOError:
            break
        
        if "**** CASE FAILED ****" in string: # parse out number of stations
            execSuccess = False
            break
        if "OVERALL PERFORMANCE" in string: # parse out number of stations
            break

    file.seek(0)
    return execSuccess

def getDimensions(file,log):
    noStations = 0 
    noSpanPts = 0 
    dataSet = False
  
    if log:
        print("attempting to retrieve dimensions. opening file",file)
  
    while not dataSet:      
        try: # read until error occurs
            string = file.readline().rstrip() # to avoid breaking on an empty line
##            print (string)
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

def saveDesign():
    
    try:
        file = open(path+inpFileName,mode='r')
    except IOError:
        print("Failed to open file"+path+inpFileName)

# load file and update data with current variable values
    newFile = [] ; counter = 0 
    while 1: 
        string = file.readline().rstrip() 
        if string == "":
            break
        else: 
            newFile.append(updateData(string)) # replace data in string with current variable values
            counter = counter + 1
    file.close()

# save the newFile
    try:
        file = open(path+inpFileName,mode='w')
    except IOError:
        print("Failed to open file: "+path+inpFileName)
    for i in range(len(newFile)):
        file.write(newFile[i]+"\n")
    file.close()

def updateData(string): 
    if string == "":
        return string
    elif string[0] == "#": 
        return string # line is a comment
    else:
        counter = 0 
        for name,variable in inpData:
            if name in string:
                element = list(inpData[counter]) # cast tuple to a list
                subString = " ".join(element[1])           
                string = name +" "+ subString
                break
            counter = counter + 1

    return string

def cleanGraph():
    plt.clf()
    main(True)
    return

def LoadOrigInpFile():
    copyfile(inpFileNameOriginal,inpFileName)
    updateModel
    python = sys.executable
    os.execl(python, python, * sys.argv)    
    return

def quitApp():
    mainWindow.destroy()
   
def readInputData(inpFileName):    
    loadData(inpFileName)

def getInpData(string,i):
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            break
        counter = counter + 1
        
    returnString = element[1]
            
    return returnString[i]

def setRobInpData(string,val,i): # robust variant of setInpData - allows number of values in tuple be increased (needs to be filled from a loop call)
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            if i >= len(element[1]): 
                element[1].append(str(round(val,5))) 
            else:
                element[1][i] = str(round(val,5)) 
            inpData[counter] = tuple(element)
            break
        counter = counter + 1
                    
    return 

def setInpData(string,val,i):
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            element[1][i] = val 
            inpData[counter] = tuple(element)
            break
        counter = counter + 1
                    
    return 

def getInpDataVal(string):
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            break
        counter = counter + 1
        
    returnString = element[1]
    
    if "." in returnString[0]: 
        return float(returnString[0])
    else: 
        return int(returnString[0])
        
        
def setInpDataVal(string,val):
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            element[1] = [str(val)] 
            inpData[counter] = tuple(element)
            break
        counter = counter + 1
                    
    return 


def getInpDataVec(string):
    counter = 0
    for name,variable in inpData:
        if name in string:
            element = list(inpData[counter]) # cast tuple to a list
            break
        counter = counter + 1
        
    stringVec = element[1]; inpVec = [0.0]*len(stringVec)
    i = 0 
    for vals in stringVec: 
        inpVec[i] = float(vals)
        i = i + 1 

    return inpVec

def getDataLabels(nData): 
    
    getDataLabels = [' ']*nData
    ptr = 0
    if nData == getNoXInFullGrid(): # 2 at inlet + interstage + 2 at outlet (canonical distribution of spline points)
        getDataLabels[ptr] = "Compressor inlet (first grid point)"
        if getInpDataVal("fl_igv:") == 1: 
            getDataLabels[ptr+1] = "IGV entry"
            getDataLabels[ptr+2] = "IGV exit - R1 entry"
            ptr = ptr + 2
        else: 
            getDataLabels[ptr+1] = "R1 entry"
            ptr = ptr + 1
            
        for i in range(getInpDataVal("n_stages:")):
            getDataLabels[ptr+1] = "R"+str(i+1)+" exit - S"+str(i+1)+" entry"
            if i+1 == getInpDataVal("n_stages:"): 
                getDataLabels[ptr+2] = "S"+str(i+1)+" exit"
            else:
                getDataLabels[ptr+2] = "S"+str(i+1)+" exit - R"+str(i+2)+" entry"                
            ptr = ptr + 2 
            
        getDataLabels[nData-1] = "Compressor exit (last grid point)"
    else: 
        for i in range(nData):
            if i == 0:
                getDataLabels[i] = "Point "+str(i+1)+" (closest to inlet)"
            else: 
                getDataLabels[i] = "Point "+str(i+1)
                
    return getDataLabels

def getFullXBezier(xmin,xmax,xHub,xTip):
    
    xHub_Bezier = []
    xHub_Bezier.append(xmin) # first point in inlet duct 
    xHub_Bezier.append(xHub[0]) # leading edge of first blade
    for i in range(1,len(xHub)-1,2):
        xHub_Bezier.append((xHub[i]+xHub[i+1])/2.0)
    xHub_Bezier.append(xHub[-1]) # trailing edge of last blade
    xHub_Bezier.append(xmax) # last point in outlet duct


    xTip_Bezier = []
    xTip_Bezier.append(xmin) # first point in inlet duct 
    xTip_Bezier.append(xTip[0]) # leading edge of first blade
    for i in range(1,len(xTip)-1,2):
        xTip_Bezier.append((xTip[i]+xTip[i+1])/2.0)
    xTip_Bezier.append(xTip[-1]) # trailing edge of last blade
    xTip_Bezier.append(xmax) # last point in outlet duct
        
    return (xHub_Bezier,xTip_Bezier)

def defineBezierWidget(sVecHubwWB,eVecHubwWB,sVecTipwWB,eVecTipwWB):

    dataLabels = getDataLabels(getInpDataVal("n_bezier:"))

    # Bezier data -----------------
    for i in range(getInpDataVal("n_bezier:")):
        Label(masterBez, text=dataLabels[i]).grid(row=i+1,column=0)
            
        sVecHubwWB.append(StringVar(masterBez)) ; sVecTipwWB.append(StringVar(masterBez)) # invoke tkinter variable constructors (StringVar) and append
        sVecHubwWB[i].set(getInpData("rbezHub:",i)) # set start variables with values contained in bezHub
        sVecTipwWB[i].set(getInpData("rbezShroud:",i)) # set start variables with values contained bezShroud
        eVecHubwWB.append(Entry(masterBez,textvariable=sVecHubwWB[i])) # use variable to set values in entries sVecHubwWB[i])
        eVecTipwWB.append(Entry(masterBez,textvariable=sVecTipwWB[i])) # use variable to set values in entries
        eVecHubwWB[i].grid(row=i+1, column=1) # arrange placement
        eVecTipwWB[i].grid(row=i+1, column=2) # arrange placement
    
    i = 0 ; val = 0
    iv = IntVar() 
    iv.set(0) # initializing the choice to the initial value of stationVariables, happens to be Vtot which happens to be first element in MODES
    Label(masterBez, text="Hub").grid(row=0, column=1)
    Label(masterBez, text="Shroud").grid(row=0, column=2)

    return

### Definition for trackGUI.py"

def postProcVarSel():
    global stationVariables
    global plotSelection
    plotSelection = str(valueRight.get())
    stationVariables= str(valueLeft.get())
    updateRunOptionsFile(path+file_name,stationVariables,plotSelection)
    file_name_2 = "dt3a.out"

    from reDraw import getRunOptions, rePlotResults
    g=plt.figure(2)
    plt.clf()
    canvas = FigureCanvasTkAgg(g, master=mainWindow)
    canvas.draw()
    canvas.get_tk_widget().grid(row=2, column=1, rowspan=4,sticky='NEWS')
    
    stationVariables, plotSelection = getRunOptions("runOptions.txt")    
    rePlotResults(stationVariables,path+file_name_2)
##    os.system("python reDraw.py") # runs reDraw.py as script
    
##def lossesSel():
##  global lossVariables
##  global plotSelection 
##  lossVariables = LOSSES[int(w.get())]
##  plotSelection = "TrailingEdges" # losses are always computed on trailing edges  
##  updateRunOptionsFile(path+file_name,lossVariables,plotSelection)
##  os.system("python reDraw.py") # runs reDraw.py as script
##

def updateRunOptionsFile(fileName,stationVariables,plotSelection):
  try:
    file = open(fileName,mode='w')
  except OSError: 
    print("failed to open file "+fileName)
    exit(0)
  
  file.write("stationVariableSelection: "+stationVariables+"\n") 
  file.write("plotSelection: "+plotSelection+"\n") 
  file.write("displayMassAveragedValue: "+str(displayMassAvgVal)+"\n") 
  file.close()
  
  return



######## End Def for trackGUI.py###

#
# retrieve data from external file

readInputData(inpFileName) # load values into inpData

rowNames = ['IGV','R1', 'S1']


mainWindow = Tk()
mainWindow.title("Interface to SC90C")
mainWindow.columnconfigure(0, weight=1)
mainWindow.columnconfigure(1, weight=3)


#getting screen width and height of display
width= mainWindow.winfo_screenwidth()
height= mainWindow.winfo_screenheight()
#setting tkinter window size
mainWindow.geometry("%dx%d" % (width, height))

# fonts for all widgets
mainWindow.option_add("*Font", "courier 8")
# font to use for label widgets
mainWindow.option_add("*Label.Font", "helvetica 8")
##mainWindow.columnconfigure(2,weight=1)



###----------------------------------------------------------------------------------------------------------###
###                                                Left Column                                               ###
###----------------------------------------------------------------------------------------------------------###

###---------- Hub / Shroud Bezier Curves sub-window ----------###

masterBez=LabelFrame(mainWindow, text="Hub/shroud bezier design def.", padx=15, pady=5)
masterBez.grid(row=0, column=0, sticky='NEWS')

defineBezierWidget(sVecHubwWB,eVecHubwWB,sVecTipwWB,eVecTipwWB) # work With Bezier data (wWB)


###-----------------------------------------------------------###

###---------- Swirl angle sub-window ----------###
masterSwirl=LabelFrame(mainWindow, text="Swirl distribution at the stator trailing edge", padx=15, pady=5)
masterSwirl.grid(row=1, column=0, sticky='NEWS')

sVecHubwWS = [] ; eVecHubwWS = [] 
sVecMidwWS = [] ; eVecMidwWS = [] 
sVecTipwWS = [] ; eVecTipwWS = []

# Stage swirl data -----------------
Label(masterSwirl, text="Hub").grid(row=0, column=1, columnspan=1)
Label(masterSwirl, text="Mid").grid(row=0, column=2, columnspan=1)
Label(masterSwirl, text="Tip").grid(row=0, column=3, columnspan=1)
if getInpDataVal("fl_igv:") == 1:
    adder = 0
else:
    adder = 1

for i in range(getInpDataVal("fl_igv:") + getInpDataVal("n_stages:")):
    if i == 0 and getInpDataVal("fl_igv:") == 1:
        Label(masterSwirl, text="IGV exit").grid(row=i+1, column=0)
    else:
        Label(masterSwirl, text="Stator " + str(i + adder) + " exit").grid(row=i + 1, column=0)
    sVecHubwWS.append(StringVar(masterSwirl))
    sVecMidwWS.append(StringVar(masterSwirl))
    sVecTipwWS.append(StringVar(masterSwirl))  # invoke tkinter variable constructors (StringVar) and append

sVecHubwWS[0].set(getInpData("alpha_TE_IGV:",0))
sVecMidwWS[0].set(getInpData("alpha_TE_IGV:",1))
sVecTipwWS[0].set(getInpData("alpha_TE_IGV:",2)) # IGV (stage 0)
sVecHubwWS[1].set(getInpData("alpha_TE_S1:", 0))
sVecMidwWS[1].set(getInpData("alpha_TE_S1:", 1))
sVecTipwWS[1].set(getInpData("alpha_TE_S1:", 2))  # Stator 1

for i in range(getInpDataVal("fl_igv:") + getInpDataVal("n_stages:")):
    eVecHubwWS.append(Entry(masterSwirl, textvariable=sVecHubwWS[i]))  # use variable to set values in entries
    eVecMidwWS.append(Entry(masterSwirl, textvariable=sVecMidwWS[i]))  # use variable to set values in entries
    eVecTipwWS.append(Entry(masterSwirl, textvariable=sVecTipwWS[i]))  # use variable to set values in entries
    eVecHubwWS[i].grid(row=i + 1, column=1)  # arrange placement
    eVecMidwWS[i].grid(row=i + 1, column=2)  # arrange placement
    eVecTipwWS[i].grid(row=i + 1, column=3)  # arrange placement

###--------------------------------------------###

###---------- Work stage data sub-window -----------###
masterWork=LabelFrame(mainWindow, width= 2, text="Radial work distribution", padx=15, pady=5)
masterWork.grid(row=2, column=0, sticky='NEWS')
sVecHubwWW = [] ; eVecHubwWW = [] 
sVecMidwWW = [] ; eVecMidwWW = [] 
sVecTipwWW = [] ; eVecTipwWW = []

# Stage work data -----------------
Label(masterWork,text="Hub").grid(row=0,column=1,columnspan=1) # top label
Label(masterWork,text="Mid").grid(row=0,column=2,columnspan=1) # top label
Label(masterWork,text="Tip").grid(row=0,column=3,columnspan=1) # top label
for i in range(getInpDataVal("n_stages:")):
    Label(masterWork, text="Rotor "+str(i+1)).grid(row=i+1,column=0)

    sVecHubwWW.append(StringVar(masterWork)) ; sVecMidwWW.append(StringVar(masterWork)) ; sVecTipwWW.append(StringVar(masterWork)) # invoke tkinter variable constructors (StringVar) and append

sVecHubwWW[0].set(getInpData("workFact_R1:",0)) ; sVecMidwWW[0].set(getInpData("workFact_R1:",1)) ; sVecTipwWW[0].set(getInpData("workFact_R1:",2)) # Rotor 1
##
for i in range(getInpDataVal("n_stages:")):
    eVecHubwWW.append(Entry(masterWork,textvariable=sVecHubwWW[i])) # use variable to set values in entries
    eVecMidwWW.append(Entry(masterWork,textvariable=sVecMidwWW[i])) # use variable to set values in entries
    eVecTipwWW.append(Entry(masterWork,textvariable=sVecTipwWW[i])) # use variable to set values in entries
    eVecHubwWW[i].grid(row=i+1, column=1) # arrange placement
    eVecMidwWW[i].grid(row=i+1, column=2) # arrange placement
    eVecTipwWW[i].grid(row=i+1, column=3) # arrange placement



###-------------------------------------------------###

###---------- Blade chord design sub-window ----------###
masterBlade=LabelFrame(mainWindow, text="Blade row design data", padx=15, pady=5)
masterBlade.grid(row=3, column=0, sticky='NEWS')

sVecwWC = [] ; eVecwWC =  [] 
CCol = 0 # Start of CCol 
sVecHubwWHS = [] ; eVecHubwWHS = [] # work With Hub and Shroud (wWHS)
sVecTipwWHS = [] ; eVecTipwWHS = [] 
# spacing
sVecwWS = [] ; eVecwWS = [] # work with Spacings (wWS)
# blades
BCol = 5 # chord stuff spans 4 columns => blade stuff start at column 4
sVecwWB = [] ; eVecwWB =  [] # work With Blades (wWB)


Label(masterBlade,text="Chord mid facts.").grid(row=0,column=CCol+1,columnspan=1)
Label(masterBlade,text="Chord hub facts.").grid(row=0,column=CCol+2,columnspan=1)
Label(masterBlade,text="Chord tip facts.").grid(row=0,column=CCol+3,columnspan=1)            
Label(masterBlade,text=  "Spacing facts.").grid(row=0,column=CCol+4,columnspan=1)

## Blade chord and number data
# chord stuff
for i in range(len(rowNames)):
    Label(masterBlade,text=rowNames[i]).grid(row=i+1,column=CCol+0)
    sVecwWC.append(StringVar()) # invoke tkinter variable constructor (StringVar) and append
    sVecwWC[i].set(getInpData("stretch:",i)) # set start variables with values contained in stretch
    eVecwWC.append(Entry(masterBlade,textvariable=sVecwWC[i])) # use variable to set values in entries
    eVecwWC[i].grid(row=i+1, column=CCol+1) # arrange placement
    sVecHubwWHS.append(StringVar()) ; sVecTipwWHS.append(StringVar()) ; sVecwWS.append(StringVar()) # invoke tkinter variable constructors (StringVar) and append
    sVecHubwWHS[i].set(getInpData("hubStretch:",i)) # set start variables with values contained in hubStretch/tipStretch
    sVecTipwWHS[i].set(getInpData("tipStretch:",i)) # set start variables with values contained in hubStretch/tipStretch
    sVecwWS[i].set(getInpData("spacing:",i)) # transfer file input for spacing: variable to sVecwWS
    eVecHubwWHS.append(Entry(masterBlade,textvariable=sVecHubwWHS[i])) # use variable to set values in entries
    eVecTipwWHS.append(Entry(masterBlade,textvariable=sVecTipwWHS[i])) # use variable to set values in entries
    eVecwWS.append(Entry(masterBlade,textvariable=sVecwWS[i])) # use variable to set values in entries
    eVecHubwWHS[i].grid(row=i+1, column=CCol+2) # arrange placement
    eVecTipwWHS[i].grid(row=i+1, column=CCol+3) # arrange placement
    eVecwWS[i].grid(row=i+1, column=CCol+4) # arrange placement

# blade stuff
Label(masterBlade,text="No Blades").grid(row=0,column=BCol+1)
for i in range(len(rowNames)):
    Label(masterBlade,text=" "+rowNames[i]).grid(row=i+1,column=BCol+0)
    sVecwWB.append(StringVar()) # invoke tkinter variable constructor (StringVar) and append
    sVecwWB[i].set(getInpData("blades:",i)) # set start variables with values contained in blades
    eVecwWB.append(Entry(masterBlade,textvariable=sVecwWB[i])) # use variable to set values in entries
    eVecwWB[i].grid(row=i+1, column=BCol+1) # arrange placement

###---------------------------------------------------###

###---------- Buttons that run SC90c ----------###
Buttons=LabelFrame(mainWindow, text="Running the case in SC90c", padx=15, pady=5)
Buttons.grid(row=4, column=0, sticky='NEWS')

Button(Buttons, text='Run', command=updateModel).grid(row=4, column=0, sticky='NEWS', padx=40, pady=4)
Button(Buttons, text='Save',   command=saveDesign).grid(row=4, column=1, sticky='NEWS',padx=40, pady=4)
Button(Buttons, text='Clean Graph',   command=cleanGraph).grid(row=4, column=2, sticky='NEWS',padx=40, pady=4)
Button(Buttons, text='Return to Original Input File',   command=LoadOrigInpFile).grid(row=4, column=3, sticky='NEWS',padx=40, pady=4)
Button(Buttons, text='Quit',   command=quitApp).grid(row=4, column=4, sticky='NEWS',padx=40, pady=4)

###---------------------------------------------###

###---------- Figure of the Compressor ----------###
CompFigure=LabelFrame(mainWindow, text="Computed Compressor Geometry", font = 'helvetica 10', padx=100, pady=0)
CompFigure.grid(row=5, column=0, sticky='NEWS')


###----------------------------------------------###

###----------------------------------------------------------------------------------------------------------###
###                                                Second Column                                             ###
###----------------------------------------------------------------------------------------------------------###

masterPost=LabelFrame(mainWindow, text="Post-processing", font = 'helvetica 10')
masterPost.grid(column=1, row=0, sticky='NEWS')

masterVariables=LabelFrame(masterPost, text="Plot Results")
masterVariables.grid(column=0, row=0, sticky='NEWS')

plots = LabelFrame(mainWindow, text="Plots", font = 'helvetica 10')
plots.grid(column=1, row=1, rowspan = 6,sticky='NEWS')


#----------------------------------------------------------------------
# Load variable names

valueLeft = StringVar()
boxLeft = ttk.Combobox(masterVariables, textvariable=valueLeft, state='readonly')
boxLeft['values'] = VARIABLES
boxLeft.current(0)
boxLeft.grid(row=0, column = 0, sticky="NEWS",padx=15)

valueRight = StringVar()
boxRight = ttk.Combobox(masterVariables, textvariable=valueRight, state='readonly')
boxRight['values'] = TYPES
boxRight.current(0)
boxRight.grid(row=0, column = 1, sticky="NEWS")

Button(masterVariables, command=postProcVarSel, text='Post-process').grid(row=0, column=2, sticky='NEWS', padx = 20)











#app=FullScreenApp(mainWindow)
mainWindow.mainloop( )


    

