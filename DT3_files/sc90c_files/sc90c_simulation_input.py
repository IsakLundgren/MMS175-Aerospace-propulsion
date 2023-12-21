import scipy.special
import scipy.optimize
import numpy
import math
import matplotlib.pyplot as plt
import sys

plotRealCompressor = False

class sc90cInput:
    def __init__(self,path,name,n_stages,fl_igv,n_span,n_duct_inlet,\
                 n_duct_outlet,n_blade_intermediate_stations,n_bezier):

# file data
        self.path = path 
        self.name = name

# basic sizing parameters
        self.n_stages = n_stages # number of stages in compressor
        self.n_span = n_span # number of spans (stream lines) in model. Must be odd number
        self.fl_igv = fl_igv # flag for IGV
        self.n_rows = 2*n_stages + fl_igv # number of blade rows
        
# stations parameters
        self.n_duct_inlet = n_duct_inlet
        self.n_duct_outlet = n_duct_outlet
        self.n_blade_intermediate_stations = n_blade_intermediate_stations
        
# additional sizing variables that are kept private for simplicity

        self.n_stations = self.n_duct_inlet + self.n_rows*2 + \
                  (self.n_rows-1)*self.n_blade_intermediate_stations + \
                     self.n_duct_outlet

# further input 
        self.Rgas = 0 # Rgas

# performance data
        self.N_rpm = 0.0    # Rotational speed [rpm]  
        self.m     = 0.0    # Mass flow [kg/s] 
        self.T_in  = 0.0    # Inlet stagnation temperature [K]
        p_in       = 0.0    # Inlet stagnation pressure [Pa]   
        
# optimization variables
        self.deltaHub = 0.0
        self.prat = 0.0
        self.reDistFact = 0.0

# pressure distribution
        self.PR_R1 = [0.0]*n_stages

# swirl distribution
        self.alpha_TE_IGV = [0.0]*n_stages
        self.alpha_TE_S1 = [0.0]*n_stages

# blade distribution
        self.noBlades = [0.0]*self.n_rows

# bezier control points
        self.n_bezier_pts = n_bezier # number of points for bezier curve
        self.x_tip_bc = [0.0]*self.n_bezier_pts ; self.x_hub_bc = [0.0]*self.n_bezier_pts 
        self.r_tip_bc = [0.0]*self.n_bezier_pts ; self.r_hub_bc = [0.0]*self.n_bezier_pts
        
# duct inlet x-data
        self.x_tip_duct_inlet = [0.0]*self.n_duct_inlet
        self.x_hub_duct_inlet = [0.0]*self.n_duct_inlet
        
# blade positions
        self.x_tip_blades = [0.0]*self.n_rows*2
        self.x_hub_blades = [0.0]*self.n_rows*2

# duct outlet x-data
        self.x_tip_duct_outlet = [0.0]*self.n_duct_outlet
        self.x_hub_duct_outlet = [0.0]*self.n_duct_outlet
        
# data holders for bezier iteration
        self.Pvec = [0.0]*self.n_bezier_pts
        self.xp = 0.0

# X and R arrays
        self.X = numpy.zeros(shape=(self.n_stations,self.n_span)) # each row holds a station
        self.R = numpy.zeros(shape=(self.n_stations,self.n_span)) # each row holds a station

# gas exit angle, degrees, or, for a rotor with NTYPE = 27, stage pressure ratio
# outlet flow angles. 
        self.RELA_bc = numpy.zeros(shape=(3,3)) # 3 rows. 3 = #(HUB;MID,TIP)
        self.RELA= numpy.zeros(shape=(self.n_rows,self.n_span)) 
        self.chord_ax = numpy.zeros(shape=(self.n_rows,self.n_span)) 


# data holders for bezier iteration
        self.x_hub = [0.0]*(self.n_duct_inlet + self.n_rows*2 + self.n_duct_outlet)
        self.x_tip = [0.0]*(self.n_duct_inlet + self.n_rows*2 + self.n_duct_outlet)


# station type indicator
        self.ISIG = [0]*self.n_stations
        self.ISTAGE = [0]*self.n_stations
        self.WALL1 = [0]*self.n_stations
        self.WALL2 = [0]*self.n_stations
        
        self.plotExists = False
        
        self.setStationIndicators()        
        
#-------------------------------------------------------------------------
#
#  Methods below
#        
# parse single line of "Axial station" record    
    def setPerformance(self,N,m,T_in,p_in):
        self.N_rpm = N; self.m = m; self.T_in = T_in; self.p_in = p_in

# parse single line of "Axial station" record            
    def setOptVars(self,prat,reDistfact):
        self.prat = prat
        self.reDistfact = reDistfact

# set pressure ratios        
    def setPressRatios(self,prat,reDistFact,rwf1):
        
        self.PR_R1 = [rwf1[0]*prat*reDistFact, rwf1[1]*prat*reDistFact  , rwf1[2]*prat*reDistFact] # [hub mid tip]
        
# set swirl distribution
    def setSwirl(self,alpha_TE_IGV,alpha_TE_S1):
        self.alpha_TE_S1  = alpha_TE_S1 # [hub mid tip]
        self.alpha_TE_IGV = alpha_TE_IGV  # [hub mid tip]

# set number of blades
    def setNoBlades(self,blades): 
        self.noBlades = blades

# define bezier curve points    
    def set_Bezier_x_Points(self,xH,xT):
        self.x_hub_bc = [xH[0]-0.5*(xH[1]-0), xH[2], xH[4], xH[5]+0.5*(xH[5]-xH[4])]
        self.x_tip_bc = [xT[0]-0.5*(xT[1]-0), xT[2], xT[4], xT[5]+0.5*(xT[5]-xT[4])]

# define bezier curve points    
    def set_Bezier_r_CurvePoints(self,bezHub,bezShr):
        self.r_hub_bc = bezHub
        self.r_tip_bc = bezShr
    
    def setDuctInlet(self,n_duct_inlet,xInlet):
        self.n_duct_inlet = n_duct_inlet
        self.x_tip_duct_inlet = xInlet 
        self.x_hub_duct_inlet = xInlet 

    def setDuctBlades(self,xHub,xTip):
        self.x_hub_blades = xHub
        self.x_tip_blades = xTip

    def setDuctOutlet(self,n_duct_outlet,xOutletHubDuct,xOutletTipDuct):
        self.n_duct_outlet = n_duct_outlet
        self.x_hub_duct_outlet = xOutletHubDuct
        self.x_tip_duct_outlet = xOutletTipDuct


    # interpolation including bezier iteration   
    def get_y_bc(self,xp,Px,Py,lastVal=5.0):
        
        self.Pvec = Px ; self.xp = xp
        
        firstVal = -1.0 #Px[0]
        
        if self.residual(firstVal)*self.residual(lastVal) >= 0.0:
            (firstVal,lastVal) = self.recoverBound(xp,Px[-1]-Px[0])
            if self.residual(firstVal)*self.residual(lastVal) >= 0.0:
                (firstVal,lastVal) = self.recoverBound(xp,Px[-1]-Px[0])
                self.assessProblem(xp,Px,Py)
                firstY = self.residual(firstVal)
                lastY  = self.residual(lastVal)
                print("This is where you will have a problem")
        
        x = scipy.optimize.brentq(self.residual,firstVal,lastVal)
        self.Pvec = Py
        get_y_bc = self.Bezier(x) # evaluate for converged x

        return get_y_bc
    
    def recoverBound(self,xp,dRange): 
        
        firstVal = xp - 0.01*dRange
        lastVal  = xp + 0.01*dRange
#
        for i in range(100): 
            if self.residual(firstVal)*self.residual(lastVal) < 0.0:
                return (firstVal,lastVal)            
            firstVal = firstVal - 0.01*dRange
            lastVal  = lastVal  + 0.01*dRange
#            f1 = self.residual(firstVal) ; f2 = self.residual(lastVal) 
                    
        return (firstVal,lastVal)            
    
    
    def assessProblem(self,xp,Px,Py):
        
        xVals = numpy.linspace(Px[0],Px[3],50).tolist()
        for i in range(len(xVals)):
            print(str(xVals[i])+" "+str(self.residual(xVals[i])))
    
    # residual
    def residual(self,t):
        
        residual = self.xp - self.Bezier(t)
        
        return residual
    
    def Bezier(self,t): 
        n = len(self.Pvec)
        
        Bezier = 0.0
        
        for i in range(n):
            val = float(scipy.special.binom(n-1,i))
            Bezier = Bezier + val*(1.0-t)**(n-1-i) * (t**i) * self.Pvec[i]
                
        return Bezier

#
# calculations
#
    def calculation(self):        
# merge lists
        self.x_hub = self.x_hub_duct_inlet + self.x_hub_blades + self.x_hub_duct_outlet # concatenate lists
        self.x_tip = self.x_tip_duct_inlet + self.x_tip_blades + self.x_tip_duct_outlet # concatenate lists

# load RELA_bc matrix
        self.RELA_bc[0] = self.alpha_TE_IGV
        self.RELA_bc[1] = self.PR_R1
        self.RELA_bc[2] = self.alpha_TE_S1
        
# x_hub and r_hub does not contain intermediate "computational" duct type stations
        r_hub = [] ; r_tip = [] 
        for i in range(len(self.x_hub)): 
            r_hub.append(self.get_y_bc(self.x_hub[i], self.x_hub_bc, self.r_hub_bc)) # interpolate and append to list
            r_tip.append(self.get_y_bc(self.x_tip[i], self.x_tip_bc, self.r_tip_bc)) # interpolate and append to list
            
        
# create X and R station representation for inlet stations
        stationCounter = 0 
        for i in range(len(self.x_hub_duct_inlet)): 
            self.X[stationCounter] = numpy.linspace(self.x_hub[i],self.x_tip[i], self.n_span)
            self.R[stationCounter] = numpy.linspace(     r_hub[i],     r_tip[i], self.n_span)
            stationCounter = stationCounter + 1


# create X and R station representation for blade stations
        for i_row in range(self.n_rows):
            LE = (2+self.n_blade_intermediate_stations)*(i_row+1) - (2+self.n_blade_intermediate_stations) + self.n_duct_inlet # leading edge pointer
            TE = (2+self.n_blade_intermediate_stations)*(i_row+1) - (1+self.n_blade_intermediate_stations) + self.n_duct_inlet # trailing edge pointer
           
            LE_data = 2*(i_row+1) - 2 + self.n_duct_inlet # leading edge pointer
            TE_data = 2*(i_row+1) - 1 + self.n_duct_inlet # trailing edge pointer 

            self.X[LE] = numpy.linspace(self.x_hub[LE_data], self.x_tip[LE_data], self.n_span) # leading edge
            self.R[LE] = numpy.linspace(     r_hub[LE_data],      r_tip[LE_data], self.n_span)

            self.X[TE] = numpy.linspace(self.x_hub[TE_data], self.x_tip[TE_data], self.n_span) # trailing edge
            self.R[TE] = numpy.linspace(     r_hub[TE_data],      r_tip[TE_data], self.n_span)
            
            r_RELA_bc = numpy.linspace(self.R[TE][0],self.R[TE][self.n_span-1],3)

            for j in range(self.n_span):
                self.RELA[i_row][j] = self.get_y_bc(self.R[TE][j],r_RELA_bc, self.RELA_bc[i_row])
                self.chord_ax[i_row][j] = math.sqrt( (self.X[TE][j] - self.X[LE][j])**2 + (self.R[TE][j] - self.R[LE][j])**2 )
        
# interblade duct stations 
        for i_row in range(self.n_stages*2+self.fl_igv-1):
            if self.n_blade_intermediate_stations ==0:
                break
            TE   = (2+self.n_blade_intermediate_stations)*(i_row+1) - (1+self.n_blade_intermediate_stations) + self.n_duct_inlet  # Trailing edge
            duct = (2+self.n_blade_intermediate_stations)*(i_row+1) - (0+self.n_blade_intermediate_stations) + self.n_duct_inlet  # Intermediate station
            LE   = (2+self.n_blade_intermediate_stations)*(i_row+1) - (-1+self.n_blade_intermediate_stations) + self.n_duct_inlet # Next leading edge
            
            for j in range(self.n_rows):
                self.X[duct][j] = (self.X[TE][j] + self.X[LE][j])/2.0

            self.R[duct][1]             = self.get_y_bc(self.X[           duct][1],self.x_hub_bc, self.r_hub_bc)
            self.R[duct][self.n_span-1] = self.get_y_bc(self.X[duct,self.n_span-1],self.x_tip_bc, self.r_tip_bc)
            self.R[duct] = numpy.linspace(self.R[duct][1],self.R[duct][self.n_span-1],self.n_span)


# outlet ducts 
        j = 0 
        for i in range(self.n_stations-self.n_duct_outlet,self.n_stations):
            duct = (self.n_duct_inlet + self.n_rows*2 + self.n_duct_outlet) - self.n_duct_outlet + j
            self.X[i] = numpy.linspace(self.x_hub[duct], self.x_tip[duct],self.n_span)
            self.R[i] = numpy.linspace(     r_hub[duct],      r_tip[duct],self.n_span)
            j = j + 1
            
        
    
    def setStationIndicators(self):        
        
# inflow duct      
        counter = 0  
        for i in range(self.n_duct_inlet): 
            self.ISIG[counter] = 0 # duct stations
            self.ISTAGE[counter] = 0 
            self.WALL1[counter] = 0 ; self.WALL2[counter] = 0 
            counter = counter + 1

# IGV
        if self.fl_igv == 1: 
            self.ISIG[counter] = 1 # IGV leading edge 
            self.ISTAGE[counter] = 0 
            self.WALL1[counter] = 0 ; self.WALL2[counter] = 0 
            counter = counter + 1            
            self.ISIG[counter] = 3 # IGV trailing edge
            self.ISTAGE[counter] = 0 
            self.WALL1[counter] = 0 ; self.WALL2[counter] = 0 
            counter = counter + 1
# intermediate duct stations
            for j in range(self.n_blade_intermediate_stations):
                self.ISIG[counter] = 0  
                self.ISTAGE[counter] = 0 
                self.WALL1[counter] = 1 ; self.WALL2[counter] = 0                 
                counter = counter + 1            
         
        for i in range(self.n_stages):
            # rotor 
            self.ISIG[counter] = 6 # Rotor leading edge 
            self.ISTAGE[counter] = i + 1 
            self.WALL1[counter] = 0 ; self.WALL2[counter] = 1                 
            counter = counter + 1            
            self.ISIG[counter] = 8 # Rotor trailing edge 
            self.ISTAGE[counter] = i + 1 
            self.WALL1[counter] = 0 ; self.WALL2[counter] = 1                 
            counter = counter + 1            
# intermediate stations
            for j in range(self.n_blade_intermediate_stations):
                self.ISIG[counter] = 0 
                self.ISTAGE[counter] = i + 1 
                self.WALL1[counter] = 1 ; self.WALL2[counter] = 0                 
                counter = counter + 1            
            
            # stator
            self.ISIG[counter] = 1 # Stator leading edge 
            self.ISTAGE[counter] = i + 1
            self.WALL1[counter] = 1 ; self.WALL2[counter] = 0                 
            counter = counter + 1            
            self.ISIG[counter] = 3 # Stator trailing edge 
            self.ISTAGE[counter] = i + 1 
            self.WALL1[counter] = 1 ; self.WALL2[counter] = 0                             
            counter = counter + 1        
            
            if i == self.n_stages - 1: 
                break
            else: 
                # intermediate stations
                for j in range(self.n_blade_intermediate_stations):
                    self.ISIG[counter] = 0 # Intermediate duct station 
                    self.ISTAGE[counter] = i + 1 
                    self.WALL1[counter] = 1 ; self.WALL2[counter] = 0                 
                    counter = counter + 1            
            
        for i in range(self.n_duct_outlet): 
            self.ISIG[counter] = 0 # duct stations
            self.ISTAGE[counter] = 0
            self.WALL1[counter] = -1 ; self.WALL2[counter] = -1                 
            counter = counter + 1

    def writeGeoFile(self,path):
        ITYPE = 0 # [0, 3, 4, 5] Can be skipped(OMITTED) => 0
        IGRID = 0 # 0 = Uniform grid. 1 = Specified radial distribution
        IMAX = self.n_stations # Number of axial planes
        IREF = self.n_duct_inlet+2 # Plane at which the grid is fixed
        ISPLIT  = IMAX + 1 # Splitter plane (if none, use IMAX+1)
        IMAXC   = IMAX # Last plane in core (if no splitter, insert IMAX)
        JMAX    = self.n_span # Number of radial stations 
        JMS     = int(self.n_span/2) # Mean streamline (it need not be the middle streamline)
        JSPLIT  = self.n_span # Number of streamline coinciding with splitter
        DATUM   = 1 # Length in metres in unit dimension

        ICURV = [0]*self.n_stations # Use 0 for a linear calculating plane (currently only valid option)


        controlParameters = [IGRID,IMAX,IREF,ISPLIT,IMAXC,JMAX,JMS,JSPLIT,DATUM]
        
        try:
            file = open(path+self.name+'.geo','w');
        except OSError:
            print("failed to open file")
            exit(0)    
          
        file.write(str(ITYPE)+"\n") 
        file.write(self.name+"\n") 
        for i in range(len(controlParameters)):
            file.write(" "+str(controlParameters[i]))
        file.write("\n") 
        for i in range(self.n_stations):
            file.write("%14.9f"% (self.X[i][0])) # X_hub
            file.write("%14.9f"% (self.X[i][self.n_span-1])) # X_shroud
            file.write("%14.9f"% (self.R[i][0])) # R_hub
            file.write("%14.9f"% (self.R[i][self.n_span-1])) # R_shroud
            file.write("%4i"% (ICURV[i])) # ICURV
            file.write("%4i"% (self.ISIG[i])) # ISIG
            file.write("%4i"% (self.ISTAGE[i])) # 
            file.write("%4i"% (self.WALL1[i])) # 
            file.write("%4i"% (self.WALL2[i])) # 
            file.write("\n") 
        file.close()


    def writeConFile(self,path):
        IFIX    = 0 # [0 - MOVING GRID], [1 - FIXED GRID]
        IREE    = 2 # [1 - Ginder formulation of radial equilibrium], [2 - Denton formulation of radial equilibrium]
        IAWBL   = 2 # [0 - OMIT W&M annulus wall boundary layer], [1 - to include it but not override AWB], [2 - include and overwrite AWB]
        IHSPS   = 1 # [1 - for H_stat procedure in inner loop], [2 - for P_stat procedure in inner loop]
        ERRM    = 0.0001 # error tolerance
        MAXCYC  = 1000 # maximum number of outer loop
        RLXFAX  = 0.8 # RELAXATION FACTOR
        RMIX    = 0.5 # relaxation factor used for mixing terms
        IPLT    = 1 # [0 - normal output files], [1 - generates the .plt file],[2 - generates the .plt but does not close the window]
        IMISES  = 0 # [0 - no mises file], [ises.xxx], [xxx.ises]
        IHAC    = 0 # [0 - no Howell & Calvert input file generated], [1 - generates the input file]
        IPRINT  = str(self.n_stations)+'*2' #  
        JPRINT  = str(self.n_span)+'*1' #  
        NCYCPR  = [0,0,0] # Three iteration numbers at which putput is to be printed.

        try:
            file = open(path+self.name+'.con','w');
        except OSError:
            print("failed to open file")
            exit(0)    

        file.write(self.name+"\n") # TITLE
        #
        file.write("%4i"% (IFIX)) 
        file.write("%4i"% (IREE)) 
        file.write("%4i"% (IAWBL)) 
        file.write("%4i"% (IHSPS)) 
        file.write("%8.5f"% (ERRM))
        file.write("%7i"% (MAXCYC))
        file.write("%6.2f"% (RLXFAX))
        file.write("%6.2f"% (RMIX))
        file.write("%4i"% (IPLT)) 
        file.write("%4i"% (IMISES)) 
        file.write("%4i"% (IHAC)) 
        file.write("\n") 
        file.write(IPRINT+"\n")
        file.write(JPRINT+"\n")
        for i in range(len(NCYCPR)):
            file.write("%4i"% (NCYCPR[i])) 
        file.write("\n") 

        file.close()
                
    def writeDefinitionFile(self,path):

        try:
            file = open(path+'sc90.fil','w');
        except OSError:
            print("failed to open file")
            exit(0)    

        file.write(self.name+".geo\n") # TITLE
        file.write(self.name+".con\n") # TITLE
        file.write(self.name+".air\n") # TITLE
        file.write(self.name+".out\n") # TITLE
        file.write("sc90c.exe\n") # TITLE
        file.close()

    def writeAirFile(self,path):

        REVS    = self.N_rpm # [rpm] Rotational speed 
        STAGT   = self.T_in # Inlet stagnation temp (K) set to 0.0 if non-uniform
        STAGP   = self.p_in # Inlet stagnation pressure (kPa) set to 0.0 if non-uniform
        ALPHA1  = 0 # inlet swirl angle (degress, positive in direction of rotor rotation) set to 100.0 if non-uniform
        TURB1   = 0 # (5) turbulence level (%) at inlet to first row (not used currently)
        TURB2   = 0 # (5) turbulence level (%) (not used)
        CP      = 0 # specific heat (kJ/(Kg K)) or gamma. Gamma assumed if 0<CP<10.0. If CP = zero, SC90 will calculate variable specific heat for air
        Rgas    = self.Rgas # gas constant (if R = 0, value of 287.054 will be used)
        Q       = 0 # fuel air ratio 
        IBLEED  = 0 # 0 if no distributed bleed, 1 for in-bleed
        JBLEED  = 0 # 0 if no wall bleed, 1 for wall in-bleed
        GAL     = 0.00056 # Gallimore mixing parameter
        CF      = 0.007 # Skin friction coefficient 
#
        NTYPE   = 27 # for design using Wright and Miller correlations, specifiying stage p.r. for rotor and stator gas outlet angles
        IANGLE  = 2 # DEFINITION OF FLOW ANGLES (1 - atan(Vw/Vmerid) 2 - atan(Vw/Vaxial))

# Wright and Miller annulus wall boundary calculation
# If any of these quantities is not known,insert zero, and a default value
# will be used
        TH1 = 0 # hub annulus wall boundary layer momentum thickness/span at first plane
        SF1 = 0 # hub annulus wall boundary layer shape factor at first plane
        TH2 = 0 # casing annulus wall boundary layer shape factor at first plane
        SF2 = 0 # casing wall boundary layer shape factor at first plane
        
        NLOSS = self.n_span # number radial positions at which data is tabulated
        NDM = NTYPE # use Miller and Wright for each rotor trailing edge
        AWB = 0 # Annulus blockage 

        TCVYSP_stator = [0.003,0.005,0.0075]+[0.0075]*self.n_stations # Tip clearance/span (assumed to be at hub for a stator, casing for a rotor) 
        TCVYSP_rotor  = [0.003,0.005,0.0075]+[0.0075]*self.n_stations# Tip clearance/span (assumed to be at hub for a stator, casing for a rotor) 

        DINC = 100.0*numpy.ones(shape=(self.n_stations,self.n_span)) # Incidence, degrees [This may be set to 100 => Wright & Miller minimum loss incidence will be used   
        BINCSS = 100.0*numpy.ones(shape=(self.n_stations,self.n_span)) # BLADE SUCTION SURFACE INCIDENCE ANGLE AT LEADING EDGE. If not known inser 100 => calculated.
        
        TOVERC_r = numpy.linspace(0.05,0.02,self.n_span)        
        TOVERC_s = numpy.linspace(0.02,0.05,self.n_span)        
        
        XFCHRD = 10.0*numpy.ones(shape=(self.n_stations,self.n_span)) # THROAT/PITCH RATIO FOR NON-DCA BLADING. (0.5 FOR DCA) (10 OR MORE AND DCA THROAT WILL BE USED)
        RLE = 0.0005*numpy.ones(shape=(self.n_stations,self.n_span)) # RADIUS OF BLADE LEADING (AND TRAILING) EDGE


        FMTOT   = self.m # inlet mass flow 
        IQ1 = 1 # plane at which core flow function to be applied
        IQ2 = 0 # plane number at which bypass flow function to be applied

        QFUNC = 0 # TOTAL FLOW FUNCTION
        QRELAX = 0 # RELAXATION FACTOR FOR CHANGES OF INLET MASS FLOW 

#        ICORRELN = 1 # CONTROL FLAG FOR MILLER AND WRIGHT CORRELATION
#        PLOSS_MAX = 0.03
#        TLOSS_MAX = 0.05
              
        try:
            file = open(path+self.name+'.air','w');
        except OSError:
            print("failed to open file")
            exit(0)    

        file.write(self.name+"\n") # TITLE
#
        file.write("{:<11.6f}\t{:>9.6f}\t{:>9.6f}\t{:2d}\t{:3d}\t{:3d}\t{:3d}".format(REVS,STAGT,STAGP/1000.0,ALPHA1,TURB1,TURB2,CP))
        file.write("{:11.3f}{:4d}{:4d}{:4d}{:12.6f}{:12.6f}\n".format(Rgas,Q,IBLEED,JBLEED,GAL,CF))
#
#       file.write("{:<6d}{:<6d}{:<6d}{:<9.6f}{:<9.6f}\n".format(NTYPE,IANGLE,ICORRELN,PLOSS_MAX,TLOSS_MAX))
        file.write("{:<6d}{:<6d}\n".format(NTYPE,IANGLE))
#
        file.write("{:<6d}{:<6d}{:<6d}{:<6d}\n".format(TH1,SF1,TH2,SF2))

#
        i_row = 0; i_rotor = 0; i_stator = 0;
        for i in range(self.n_stations): # 0 - for duct region, 1 - Stator LE, 3 - Stator TE, 6 - Rotor LE, 8 Rotor TE
            if self.ISIG[i] == 3: # stator
#
                file.write("%4i"% (i+1)) 
                file.write("\n") 
#                
                file.write("%4i"% (NLOSS)) 
                file.write("%4i"% (NDM)) 
                file.write("%4i"% (AWB)) 
                file.write("\n") 
                
                file.write("%4i"% (self.noBlades[i_row])) 
                file.write("%9.4f"% (TCVYSP_stator[i_stator]))
                file.write("\n") 
                
                for j in range(self.n_span):                    
                    file.write("%12.7f"% (self.R[i-1][j]))
                    file.write("%12.7f"% (self.R[i][j]))
                    file.write("%8.2f"% (DINC[i_row][j]))
                    file.write("%8.2f"% (BINCSS[i_row][j]))
                    file.write("%12.7f"% (self.RELA[i_row][j]))
                    file.write("%12.7f"% (TOVERC_s[j]))
                    file.write("%12.7f"% (XFCHRD[i_row][j]))
                    file.write("%12.7f"% (RLE[i_row][j]))
                    file.write("%12.7f"% (self.chord_ax[i_row][j]))
                    file.write("\n")        
         
                i_row   = i_row + 1
                i_stator = i_stator + 1
            
            elif self.ISIG[i] == 8: 
                #
                file.write("%4i"% (i+1))
                file.write("\n") 
                #                
                file.write("%4i"% (NLOSS)) 
                file.write("%4i"% (NDM)) 
                file.write("%4i"% (AWB)) 
                file.write("\n") 
                #
                file.write("%4i"% (self.noBlades[i_row])) 
                file.write("%9.4f"% (TCVYSP_rotor[i_rotor]))
                file.write("\n") 

                for j in range(self.n_span):                    
                    file.write("%12.7f"% (self.R[i-1][j]))
                    file.write("%12.7f"% (self.R[i][j]))
                    file.write("%8.2f"% (DINC[i_row][j]))
                    file.write("%8.2f"% (BINCSS[i_row][j]))
                    file.write("%12.7f"% (self.RELA[i_row][j]))
                    file.write("%12.7f"% (TOVERC_s[j]))   # this is probably wrong
                    file.write("%12.7f"% (XFCHRD[i_row][j]))
                    file.write("%12.7f"% (RLE[i_row][j]))
                    file.write("%12.7f"% (self.chord_ax[i_row][j]))
                    file.write("\n")        

                i_row   = i_row + 1
                i_rotor = i_rotor + 1
                
        file.write("%4i"% (0))
        file.write("\n") 
        #
        file.write("%9.4f"% (FMTOT))
        file.write("%9.4f"% (IQ1)) 
        file.write("%9.4f"% (IQ2)) 
        file.write("%9.4f"% (QFUNC)) 
        file.write("%9.4f"% (QRELAX)) 
        file.write("\n") 
        #
        file.write(str(self.n_stations)+"*1.0")   # CH - 1.0 for subsonic solution at plane I

        file.close()


    def createSC90C_Input_Files(self):
        
        self.writeDefinitionFile(self.path)

        
        self.writeGeoFile(self.path)
        self.writeConFile(self.path)
        self.writeAirFile(self.path)
        
        
    def aBladeStation(self,i):
        
        aBladeStation = False    
        if self.ISIG[i] > 0: 
            aBladeStation = True
        
        return aBladeStation

    
    def getBladeNumbers(self,rowType,nbl):
        
        getBladeNumbers = []
        ptr = 0 
        if self.fl_igv == 1 and rowType=="stator":
            getBladeNumbers.append(nbl[ptr])
            ptr = ptr + 1
        
        for i in range(self.fl_igv,2*self.n_stages,2): 
            if rowType == "rotor": 
                getBladeNumbers.append(nbl[i])
            else: 
                getBladeNumbers.append(nbl[i+1])
            ptr = ptr + 1
        
        return getBladeNumbers
    
    def writeVec(self,file,str,vec): 
        
        file.write(str) 
        for i in range(len(vec)): 
            file.write("%8.4f"% (vec[i]))    
        file.write("\n")
        
        return 
    
    def writeiVec(self,file,str,vec): 
        
        file.write(str) 
        for i in range(len(vec)): 
            file.write("%5i"% (vec[i]))    
        file.write("\n")
        
        return 

    def createPlots(self,streamLines):
        import matplotlib.pyplot as plt
        plt.style.use('seaborn-paper')
        axialLength = self.X[self.n_stations-1,0] - self.X[0,0]
        xMax = self.X[self.n_stations-1,0]; xMin = self.X[0,0]

        radialHeight = self.R[0,self.n_span-1] - self.R[0,0]
        rMax = self.R[0,self.n_span-1] ; rMin = self.R[0,0]

        # plot hub and shroud line
        plt.plot(self.X[:,0],self.R[:,0],'o-',color="black",linewidth=1.0)
        plt.plot(self.X[:,self.n_span-1],self.R[:,self.n_span-1],'o-',color="black",linewidth=1.0)


        if plotRealCompressor:
            Rt = self.R[self.n_duct_inlet+2*self.fl_igv,self.n_span-1]  # Rtip for rotor leading edge
            x0 = self.X[self.n_duct_inlet,0]
            xhub, yhub, xtip, ytip = self.retrieveDataFromFile(Rt,x0,"scannedGasPath.txt")
            plt.plot(xhub,yhub,'ro-',linewidth=2.0)
            plt.plot(xtip,ytip,'ro-',linewidth=2.0)
            for i in range(len(xhub)):
                x = [xhub[i],xtip[i]] ; y = [yhub[i],ytip[i]]
                plt.plot(x,y,'ro-',linewidth=1.0)

            htr = yhub[3]/ytip[3]

        # plot stream lines
        for i in range(1,self.n_span-1):
            xvec = []
            for j in range(self.n_stations):
                xvec.append(streamLines[j].variables[0][i])
            rvec = []
            for j in range(self.n_stations):
                rvec.append(streamLines[j].variables[1][i])
            plt.plot(xvec,rvec,color="blue",linewidth=0.3)


        # plot leading edge and trailing edge
        for i in range(self.n_stations):
            if self.aBladeStation(i):
                plt.plot(self.X[i,:],self.R[i,:],'-',color="red",linewidth=1.0)

        plt.xlabel('Axial coordinate (m)')
        plt.ylabel('Radial coordinate (m)')
        plt.title('Duct geometry')
        plt.grid(True)

        #Plot bezier
        plt.plot(self.x_tip_bc, self.r_tip_bc, 's--',linewidth=0.2, label='Bezier control points - shroud')
        plt.plot(self.x_hub_bc, self.r_hub_bc, 's--',linewidth=0.2, label='Bezier control points - hub')

        plt.legend()
        if axialLength > radialHeight:
            plt.axis('equal')
            plt.axis([xMin-0.10*axialLength,xMax+0.10*axialLength,rMax-0.8*axialLength,rMax+0.40*axialLength])
##        plt.show(block=False)
        
        
    def retrieveDataFromFile(self,Rt,x0,fname): 
        tags = ["hub:","tip:","scale:","#"]
        
        try:
            file = open(self.path+fname,'r');
        except OSError:
            print("failed to open file: "+self.path+fname)
            exit(0)    
        
        complete = False
        xhub = [] ; yhub = [] ; xtip = [] ; ytip = []
        while not complete:
            try: # read until error occurs
                string = file.readline().rstrip() # to avoid breaking on an empty line
            except IOError:
                break
            
            if "hub:" in string: # parse out number of stations
                xhub, yhub = self.getData(tags,file)
            
            if "tip:" in string: # parse out number of stations
                xtip, ytip = self.getData(tags,file)
                
            if "scale:" in string: 
                string = file.readline().rstrip() # numbers are on next line
                subStrings = string.split()
                xscale = float(subStrings[0])
                yscale = float(subStrings[1])
                complete = True

        # [0] - IGV leading edge, [1] - IGV trailing edge, [2] - Rotor leading edge in data set
        compScale = 1.05 ; scanScale = 1.0 # yscale / xscale
        for i in range(len(xtip)): # make image properly scaled - not distorted
            xtip[i] = scanScale*xtip[i]
            xhub[i] = scanScale*xhub[i]
            
        for i in range(len(xtip)): # scale to real compressor size
            xtip[i] = compScale*xtip[i]
            xhub[i] = compScale*xhub[i]
            ytip[i] = compScale*ytip[i]
            yhub[i] = compScale*yhub[i]
            
        # translate in x-coordinate
        xdelta = xhub[0]-x0
        for i in range(len(xtip)): # scale to real compressor size
            xtip[i] = xtip[i] - xdelta
            xhub[i] = xhub[i] - xdelta
            
        
        return (xhub, yhub, xtip, ytip)
    
    def getData(self,tags,file): 
        
        x = [] ; y = [] 
        while 1:
            string = file.readline().rstrip() # to avoid breaking on an empty line
            if self.isATag(tags,string): 
                break
            else:
                subStrings = string.split()
                x.append(float(subStrings[0])) ; y.append(float(subStrings[1]))
                
        return (x, y)

    def isATag(self,tags,string): 
        
        isATag = False
        for str in tags: 
            if str in string:
                isATag = True
                break
            
        return isATag

# xTip is specified at both leading and trailing edge. noPts = (2*stages + IGV)*2 => 14 for VINK
    def stretchChord(self,xTip,stretch,spacing): 
        
# if all elements in spacing are one, xTip is left unaffected
        bladeCounter = 0 ; noX = len(xTip) ; deltaTot = 0.0
        for i in range(2,noX): # first two are first blade => no upstream spacing.
            if i % 2 == 0: # if i is even => leading edge. (0,2,4,6,8,... are leading edges)
                delta = (xTip[i] - xTip[i-1])*(spacing[bladeCounter]-1.0)
                for j in range(i,noX): # add delta to all downstream x-coordinates
                    xTip[j] = xTip[j] + delta
                deltaTot = deltaTot + delta
                bladeCounter = bladeCounter + 1


# if all elements in strech are one, xTip is left unaffected
        bladeCounter = 0 ; noX = len(xTip) ; deltaTot = 0.0
        for i in range(noX): 
            if i % 2 == 0: # if i is even => leading edge. (0,2,4,6,8,... are leading edges)
                delta = (xTip[i+1] - xTip[i])*(stretch[bladeCounter]-1.0)
                for j in range(i+1,noX): # add delta to all downstream x-coordinates
                    xTip[j] = xTip[j] + delta
                deltaTot = deltaTot + delta
                bladeCounter = bladeCounter + 1
        
        return (xTip, deltaTot)     
    
    def isALeadingEdge(self,stationType): 
        isALeadingEdge = False
        if stationType == 1 or stationType == 6: 
            isALeadingEdge = True
            
        return isALeadingEdge
    
    def logAspectRatios(self): 

# log stator aspect ratios
        sys.stdout.write('Stator aspect ratios are: ')
        
        ars = [] 
        for i in range(self.n_stations): 
            if self.ISIG[i] == 1: 
                
                xHubLe = self.X[i][0] ; xHubTe = self.X[i+1][0]
                xTipLe = self.X[i][self.n_span-1] ; xTipTe = self.X[i+1][self.n_span-1]                
                rHubLe = self.R[i][0] ; rHubTe = self.R[i+1][0]
                rTipLe = self.R[i][self.n_span-1] ; rTipTe = self.R[i+1][self.n_span-1]
#                
                chord  = ( (xHubTe-xHubLe)  + (xTipTe-xTipLe) )/2.0
                height = ( (rTipLe+rTipTe)/2.0 - (rHubLe+rHubTe)/2.0 )
                
                sys.stdout.write("{:6.2f}".format(height/chord))
                ars.append(height/chord)
        print("\n")

        sys.stdout.write('Rotor aspect ratios are: ')

        arr = [] 
        for i in range(self.n_stations): 
            if self.ISIG[i] == 6: 
                
                xHubLe = self.X[i][0] ; xHubTe = self.X[i+1][0]
                xTipLe = self.X[i][self.n_span-1] ; xTipTe = self.X[i+1][self.n_span-1]                
                rHubLe = self.R[i][0] ; rHubTe = self.R[i+1][0]
                rTipLe = self.R[i][self.n_span-1] ; rTipTe = self.R[i+1][self.n_span-1]
                #
                chord  = ( (xHubTe-xHubLe)  + (xTipTe-xTipLe) )/2.0
                height = ( (rTipLe+rTipTe)/2.0 - (rHubLe+rHubTe)/2.0 )
                
                sys.stdout.write("{:6.2f}".format(height/chord))

                arr.append(height/chord)

        print("\n")

        return (ars,arr)
    

