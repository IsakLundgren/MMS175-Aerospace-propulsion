import os, glob
import numpy as np
from scipy.optimize import fsolve
import datetime
import math
import shutil


#================ PYTHON SCRIP FOR saveCompressorForSmooth ================
# Derived from saveCompressorForSmooth.m (MatLab file)
# Feel free to save and use for next years students (python users) 
# Composed by:
#       Robert Ranman
#       Chalmers Univesity of Technology
#       2021-12-12
#==========================================================================

CWD = os.path.dirname(os.path.realpath(__file__))

def logVec(file, string, npts, vec):
    file.write("%s" % (string))
    for i in range(npts):
        file.write("%8.4f" % (vec[i]))
    file.write("\n")

def nchoosek(n, k):
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)

def get_B_curve(P, t):
    n = len(P)
    B = 0
    for i in range(n):
        B += nchoosek(n-1, i) * (1-t)**(n-1-i) * t**i * (P[i])
    return B

def get_y_bc(x_point, Px, Py):
    t_in = (x_point - Px[0])/(Px[-1] - Px[0]) #initial guess
    t = fsolve(lambda t: x_point - get_B_curve(Px,t), t_in, factor=0.1)
    y_point = get_B_curve(Py, t)
    return y_point


def getBezierData(n_bezier_pts, x_igv_le, x_igv_te, x_rot_le, x_stat_le, x_stat_te, r_igv_tle, r_rot_tte, r_stat_tle, r_stat_tte, r_igv_hle, r_rot_hte, r_stat_hle, r_stat_hte):

    if n_bezier_pts == 4:
        # define xpts for hub and shroud bezier curve
        axChord = x_igv_te - x_igv_le                 #
        x1 = x_igv_le - 0.5*axChord                   # one blade chord upstream of igv leading edge
        x11 = x_rot_le
        x2 = x_stat_le                                 # igv trailing edge
        x3 = x_stat_te + 0.5*(x_stat_te - x_stat_le)  # HPCgeom.Stator.xte(1) + 1.0*(HPCgeom.Stator.xte(1)-HPCgeom.Stator.xle(1)); % one stator chord downstream of S1 trailing edge
        # number of x-stations above must equal n_bezier_pts
   
        xbezHub = [x1, x11, x2, x3]
        xbezShroud = xbezHub # no sweep is assumed

        # bezier interpolation for shroud lines 
        r1 = r_igv_tle   #constant radius one blade chord upstream of IGV leading edge
        r11 = r_igv_tle
        r2 = r_rot_tte
        r3 = r_stat_tte  #  constant radius one stator chord downstream of S1 trailing edge

        # shroud interpolation to smoothen wall contour
        rbezShroud = np.zeros((n_bezier_pts, 1))
        for i in range(n_bezier_pts): 
            rbezShroud[i] = get_y_bc(xbezShroud[i],xbezShroud,[r1,r11, r2, r3])


        # bezier interpolation for hub lines 
        r1 = r_igv_hle     # constant radius one blade chord upstream of igv leading edge
        r11 = r_rot_hte
        r2 = r_rot_hte     # igv trailing edge
        r3 = r_stat_hte    # constant radius one stator chord downstream of S1 trailing edge

        # hub interpolation to smoothen wall contour
        rbezHub = np.zeros((n_bezier_pts, 1))
        for i in range(n_bezier_pts): 
            rbezHub[i] = get_y_bc(xbezHub[i],xbezHub,[r1, r11, r2, r3])

    return xbezHub, rbezHub, xbezShroud, rbezShroud

def saveCompressorForSmooth(NH, m, tin, pin, pout,x_igv_le, x_igv_te, x_rot_le, x_rot_te, x_stat_le, x_stat_te,
                            r_igv_tle, r_igv_tte, r_rot_tle, r_rot_tte, r_stat_tle, r_stat_tte, r_igv_hle, r_igv_hte,r_rot_hle, r_rot_hte, r_stat_hle, r_stat_hte):

    #============================================
    #   Saving Backup Fils for smoothInput.txt
    #============================================
    # if glob.glob(os.path.join(CWD, "sc90c_files/smoothInput.txt")):
    #     if glob.glob(os.path.join(CWD, "sc90c_files/backup")):
    #         FilePath = os.path.join(CWD, "sc90c_files\\smoothInput.txt")
    #         created = os.path.getmtime(FilePath)
    #         timestamp = datetime.datetime.fromtimestamp(created)
    #         os.rename(FilePath, os.path.join(CWD, "sc90c_files/backup/smoothInput_%s.txt" % (str(timestamp).replace(":", "-").replace(" ", "_").split(".")[0]))) 
    #         print("Saving Backup File...")
    #     else:
    #         os.mkdir(os.path.join(CWD, "sc90c_files/backup"))
    #         FilePath = os.path.join(CWD, "sc90c_files\\smoothInput.txt")
    #         created = os.path.getmtime(FilePath)
    #         timestamp = datetime.datetime.fromtimestamp(created)
    #         os.rename(FilePath, os.path.join(CWD, "sc90c_files/backup/smoothInput_%s.txt" % (str(timestamp).replace(":", "-").replace(" ", "_").split(".")[0]))) 
    #         print("Saving Backup File...")
 
    prat = pout/pin

    with open(os.path.join(CWD, "sc90c_files/smoothInput.txt"), "w") as f:
        f.write("%s %4.4f \n" % ("prat:", prat))
        f.write("%s %4.4f \n" % ("NH:", 60.0*NH))
        f.write("%s %4.4f \n" % ("mdot:", m))
        f.write("%s %4.4f \n" % ("T0:", tin))
        f.write("%s %4.4f \n" % ("P0:", pin))
        f.write("%s %4.4f \n" % ("reDistFact:", 1.0))
        f.write("%s %i \n" % ("n_stages:", 1))
        f.write("%s %i \n" % ("n_span:", 11))
        f.write("%s %i \n" % ("fl_igv:", 1))


        n_bezier_pts = 4

        f.write("%s %i \n" % ("n_bezier:", n_bezier_pts))        
        f.write("%s %i \n" % ("n_duct_inlet:", 2))
        f.write("%s %i \n" % ("n_duct_outlet:", 2))
        f.write("%s %i \n" % ("n_blade_intermediate_stations:", 0))
        logVec(f, "hubStretch: ",3 ,np.ones((3)))
        logVec(f, "tipStretch: ",3 ,np.ones((3)))
        logVec(f, "stretch: ",3 ,np.ones((3)))
        logVec(f, "spacing: ",3 ,np.ones((3)))
        logVec(f, "hubSpacing: ",3 ,np.ones((3)))
        logVec(f, "tipSpacing: ",3 ,np.ones((3)))

        f.write("%s %i %i %i \n" % ("blades: ", 70, 100, 150))

        
        xbezHub, rbezHub, xbezShroud, rbezShroud = getBezierData(n_bezier_pts, x_igv_le,x_igv_te, x_rot_le, x_stat_le,
                                                                 x_stat_te, r_igv_tle, r_rot_tte, r_stat_tle,
                                                                 r_stat_tte, r_igv_hle, r_rot_hte, r_stat_hle, r_stat_hte)

        logVec(f, "xbezHub: ", n_bezier_pts, xbezHub)
        logVec(f, "rbezHub: ", n_bezier_pts, rbezHub)
        logVec(f, "xbezShroud: ", n_bezier_pts, xbezShroud)
        logVec(f, "rbezShroud: ", n_bezier_pts, rbezShroud)

        logVec(f, "alpha_TE_S1: ", 3, [0.0,  0.0,  0.0])
        logVec(f, "alpha_TE_IGV: ", 3, [0.0, 0.0, 0.0])

        f.write("%s %4.4f %4.4f %4.4f  \n" % ("workFact_R1: ", 1, 1, 1))


        logVec(f, "xHub: ", 6, [x_igv_le, x_igv_te, x_rot_le, x_rot_te, x_stat_le, x_stat_te])
        logVec(f, "xTip: ", 6, [x_igv_le, x_igv_te,x_rot_le, x_rot_te, x_stat_le, x_stat_te])

        f.close()
        shutil.copy2(os.path.join(CWD, "sc90c_files/smoothInput.txt"), os.path.join(CWD, "sc90c_files/smoothInputOriginal.txt"))

