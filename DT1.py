import numpy as np
import matplotlib.pyplot as plt


# Given parameters
Altitude = 35000 * 0.3048  # m
Flight_Mach_number = 0.82
dT_isa = 10  # K
Net_thrust = 67.3 * 1e3  # N
BPR = 11.8  # Bypass ratio
FPR = 1.55  # Fan pressure ratio
OPR = 47  # Overall pressure ratio
HPC_pressure_ratio = 4.5
Turbine_inlet_temperature = 1680  # K
Intake_efficiency = 0.995  #
Fan_polytropic_efficiency = 0.92  #
Combustor_efficiency = 0.99  #
Combustor_pressure_loss = 0.035  #
HPT_polytropic_efficiency = 0.905  #
IPT_polytropic_efficiency = 0.91  #
LPT_polytropic_efficiency = 0.915  #
Cold_Jet_efficiency = 0.98  #
Hot_jet_efficiency = 0.99  #
Shaft_mechanical_efficiency = 0.995  #
IPC_polytropic_efficiency = 0.91  #
HPC_polytropic_efficiency = 0.915  #

# Thermodynamic ambient properties at design altitude
sourceAlt = [10000, 11000]  # m https://www.engineeringtoolbox.com/elevation-speed-sound-air-d_1534.html
sourceTemperature = [-49.9 + 273.15, -56.4 + 273.15]  # K
sourcePressure = [26.48 * 1e3, 22.68 * 1e3]  # Pa
sourceSOS = [299.5, 295.2]  # m/s
T_1 = np.interp(Altitude, sourceAlt, sourceTemperature) + dT_isa   # K
p_1 = np.interp(Altitude, sourceAlt, sourcePressure)  # Pa
gamma_a = 1.4
gamma_g = 1.333
cp_a = 1005  # J/kgK
cp_g = 1148  # # J/kgK
R_a = 287     # unit
R_g = cp_g * (gamma_g-1)/gamma_g
Speed_of_sound_1 = np.sqrt(gamma_a * R_a * T_1)  # m/s

# Chemical energy in jet-A
Q_net_JA = 43.1e6   # J/kg

# Flight velocity
C_a = Flight_Mach_number * Speed_of_sound_1  # m/s

# IPC pressure ratio
IPC_pressure_ratio = OPR / (FPR * HPC_pressure_ratio)


def massFlowToThrust(dmdt_0, coolingFraction=0.0, coolSplitFrac=0.0, hasPrinting=False, printLatex=False,
                     BPR=11.8, FPR=1.55):
    assert(1 > coolingFraction >= 0)
    assert(1 > coolSplitFrac >= 0)
    hasCooling = coolingFraction > 0

    dmdt_hot = dmdt_0 / (BPR + 1)
    dmdt_cold = BPR * dmdt_hot

    # Follow thermodynamic quantities along flow path
    # Note that the measurements take place directly after the component
    # Assume calorically perfect gas
    T0_1 = T_1 * (1 + (gamma_a - 1) / 2 * Flight_Mach_number ** 2)
    p0_1 = p_1 * (1 + (gamma_a - 1) / 2 * Flight_Mach_number ** 2) ** (gamma_a / (gamma_a - 1))

    p0_2 = p0_1 * FPR   # pressure after fan
    T0_2 = T0_1 * FPR ** ((gamma_a - 1) / (Fan_polytropic_efficiency * gamma_a))
    dWdt_fan = dmdt_0 * cp_a * (T0_2 - T0_1)

    p0_3 = p0_2 * IPC_pressure_ratio    # pressure after IPC
    T0_3 = T0_2 * IPC_pressure_ratio ** ((gamma_a - 1) / (IPC_polytropic_efficiency * gamma_a))
    dWdt_IPC = dmdt_hot * cp_a * (T0_3 - T0_2)

    p0_4 = p0_3 * HPC_pressure_ratio    # pressure after HPC
    T0_4 = T0_3 * HPC_pressure_ratio ** ((gamma_a - 1) / (HPC_polytropic_efficiency * gamma_a))
    dWdt_HPC = dmdt_hot * cp_a * (T0_4 - T0_3)

    # fuel air mass ratio
    FAR_alpha = 0.10118 + 2.00376E-05 * (700 - T0_4)
    FAR_beta = 3.7078E-03 - 5.2368E-06 * (700 - T0_4) - 5.2632E-06 * Turbine_inlet_temperature
    FAR_gamma = 8.889E-08 * (Turbine_inlet_temperature - 950)
    FAR = (FAR_alpha - np.sqrt(FAR_alpha ** 2 + FAR_beta) - FAR_gamma) / Combustor_efficiency

    dmdt_CC = (1 - coolingFraction) * dmdt_hot
    dmdt_cool = coolingFraction * dmdt_hot
    dmdt_cool_s = (1 - coolSplitFrac) * dmdt_cool
    dmdt_cool_r = coolSplitFrac * dmdt_cool
    dmdt_f = FAR * dmdt_CC
    dmdt_g_1 = dmdt_f + dmdt_CC

    T0_5 = Turbine_inlet_temperature
    p0_5 = p0_4 * (1 - Combustor_pressure_loss)

    dmdt_g_2 = dmdt_g_1 + dmdt_cool_s
    T0_5_2 = (dmdt_g_1 * cp_g * T0_5 + dmdt_cool_s * cp_a * T0_4) / (cp_g * dmdt_g_2)

    T0_5_3 = T0_5_2 - dWdt_HPC / (dmdt_g_2 * cp_g * Shaft_mechanical_efficiency)
    dmdt_g = dmdt_g_2 + dmdt_cool_r
    T0_6 = (dmdt_g_2 * cp_g * T0_5_3 + dmdt_cool_r * cp_a * T0_4) / (cp_g * dmdt_g)
    p0_6 = p0_5 * (T0_5_3 / T0_5_2) ** (gamma_g / ((gamma_g - 1) * HPT_polytropic_efficiency))

    T0_7 = T0_6 - dWdt_IPC / (dmdt_g * cp_g * Shaft_mechanical_efficiency)
    p0_7 = p0_6 * (T0_7 / T0_6) ** (gamma_g / ((gamma_g - 1) * IPT_polytropic_efficiency))

    T0_8 = T0_7 - dWdt_fan / (dmdt_g * cp_g * Shaft_mechanical_efficiency)
    p0_8 = p0_7 * (T0_8 / T0_7) ** (gamma_g / ((gamma_g - 1) * LPT_polytropic_efficiency))

    # pressure ratio critical
    CPR_hot = (1 / (1 - 1 / Hot_jet_efficiency * (gamma_g - 1) / (gamma_g + 1)) ** (gamma_g/(gamma_g-1)))
    CPR_cold = (1 / (1 - 1 / Cold_Jet_efficiency * (gamma_a - 1) / (gamma_a + 1)) ** (gamma_a/(gamma_a-1)))

    hotIsChoked = (p0_8 / p_1) > CPR_hot
    coldIsChoked = (p0_2 / p_1) > CPR_cold

    if hotIsChoked:  # p9=p9c critical pressure
        p_9 = p0_8 / CPR_hot
        T_9 = 2*T0_8/(gamma_g+1)
        C_9 = np.sqrt(gamma_g*R_g*T_9)
        rho_9 = p_9/(R_g*T_9)
        A_9 = dmdt_g / (rho_9*C_9)
        F_GH = dmdt_g * C_9 + A_9*(p_9-p_1)  # Gross hot gases thrust
        hot_pratio = CPR_hot
        # ideal velocity
        C9_ideal = np.sqrt(2 * cp_g * T0_8 * (1 - (1 / (p0_8 / p_1)) ** ((gamma_g - 1) / gamma_g)))
        C_9_eta = C9_ideal
    else:
        p_9 = p_1
        T_9 = T0_8-Hot_jet_efficiency*T0_8*(1-(1/(p0_8/p_9))**((gamma_g-1)/gamma_g))
        C_9 = np.sqrt((T0_8-T_9)*2*cp_g)
        F_GH = dmdt_g * C_9
        hot_pratio = (p0_8 / p_1)
        C_9_eta = C_9

    if coldIsChoked:
        p_10 = p0_2 / CPR_cold
        T_10 = 2 * T0_2 / (gamma_a + 1)
        # velocity = speed of sound
        C_10 = np.sqrt(gamma_a * R_a * T_10)
        rho_10 = p_10 / R_a / T_10
        A_10 = dmdt_cold / rho_10 / C_10
        F_GC = dmdt_cold * C_10 + A_10 * (p_10 - p_1)
        cold_pratio = CPR_cold
        # ideal velocity
        C_10_ideal = np.sqrt(2 * cp_a * T0_2 * (1 - (1 / (p0_2 / p_1)) ** ((gamma_a - 1) / gamma_a)))
        C_10_eta = C_10_ideal
    else:
        p_10 = p_1
        T_10 = T0_2 - Cold_Jet_efficiency * T0_2 * (1 - (1 / (p0_2 / p_10)) ** ((gamma_a - 1) / gamma_a))
        C_10 = np.sqrt((T0_2 - T_10) * 2 * cp_a)
        F_GC = dmdt_cold * C_10
        cold_pratio = (p0_2 / p_1)
        C_10_eta = C_10

    F_D = dmdt_0 * C_a

    F_net = F_GH + F_GC - F_D
    SFC = dmdt_f / F_net

    deltaW_kin = dmdt_cold * C_10_eta ** 2 / 2 + dmdt_g * C_9_eta ** 2 / 2 - dmdt_0 * C_a**2 / 2

    eta_p = F_net * C_a / deltaW_kin
    eta_th = deltaW_kin / dmdt_f / Q_net_JA
    eta_0 = eta_p * eta_th

    xi = C_10_eta / C_9_eta

    dmdt = np.array([0, dmdt_0, dmdt_hot, dmdt_hot, dmdt_CC, dmdt_CC, dmdt_g, dmdt_g, dmdt_g])
    T0 = np.array([0, T0_1, T0_2, T0_3, T0_4, T0_5, T0_6, T0_7, T0_8])
    p0 = np.array([0, p0_1, p0_2, p0_3, p0_4, p0_5, p0_6, p0_7, p0_8])
    gamma = np.array([0, gamma_a, gamma_a, gamma_a, gamma_a, gamma_g, gamma_g, gamma_g, gamma_g])
    R = np.array([0, R_a, R_a, R_a, R_a, R_g, R_g, R_g, R_g])
    cp = np.array([0, cp_a, cp_a, cp_a, cp_a, cp_g, cp_g, cp_g, cp_g])


    if hasPrinting:
        if printLatex:
            if hasPrinting:
                station_properties_table = "\\begin{table}[ht]\n"
                station_properties_table += "\\centering\n"
                station_properties_table += ("\\caption{Analysis results, cooling fraction: "
                                             + f"{coolingFraction * 100:.2g}" +
                                             "\\%.}\n")
                station_properties_table += "\\begin{subtable}{0.45\\textwidth}\n"
                station_properties_table += "\\centering\n"
                station_properties_table += "\\caption{Station Thermodynamic Properties}\n"
                station_properties_table += "\\begin{tabular}{|l|l|l|l|}\n"
                station_properties_table += "\\hline\n"
                station_properties_table += "Station & $P_0$ (kPa) & $T_0$ (K) & $\\dot{m}$ (kg/s) \\\\\n"
                station_properties_table += "\\hline\n"

                station_data = [
                    ["1", f"{p0_1 * 1e-3:.4g}", f"{T0_1:.4g}", f"{dmdt_0:.4g}"],
                    ["2", f"{p0_2 * 1e-3:.4g}", f"{T0_2:.4g}", f"{dmdt_hot:.4g}"],
                    ["3", f"{p0_3 * 1e-3:.4g}", f"{T0_3:.4g}", f"{dmdt_hot:.4g}"],
                    ["4", f"{p0_4 * 1e-3:.4g}", f"{T0_4:.4g}", f"{dmdt_hot:.4g}"],
                ]

                for data in station_data:
                    station_properties_table += " & ".join(data) + " \\\\\n"

                if hasCooling:
                    station_data_cooling = [
                        ["5.1", f"{p0_5 * 1e-3:.4g}", f"{T0_5:.4g}", f"{dmdt_g_1:.4g}"],
                        ["5.2", f"{p0_5 * 1e-3:.4g}", f"{T0_5_2:.4g}", f"{dmdt_g_2:.4g}"],
                        ["5.3", f"{p0_6 * 1e-3:.4g}", f"{T0_5_3:.4g}", f"{dmdt_g_2:.4g}"],
                    ]

                    for data in station_data_cooling:
                        station_properties_table += " & ".join(data) + " \\\\\n"
                else:
                    station_data_no_cooling = [
                        ["5.1", f"{p0_5 * 1e-3:.4g}", f"{T0_5:.4g}", f"{dmdt_g:.4g}"],
                        ["5.2", f"{p0_5 * 1e-3:.4g}", f"{T0_5_2:.4g}", f"{dmdt_g:.4g}"],
                        ["5.3", f"{p0_6 * 1e-3:.4g}", f"{T0_5_3:.4g}", f"{dmdt_g:.4g}"],
                    ]

                    for data in station_data_no_cooling:
                        station_properties_table += " & ".join(data) + " \\\\\n"

                station_properties_table += ("6 & " + f"{p0_6 * 1e-3:.4g}" + " & " + f"{T0_6:.4g}" +
                                             " & " + f"{dmdt_g:.4g}" + " \\\\\n")
                station_properties_table += ("7 & " + f"{p0_7 * 1e-3:.4g}" + " & " + f"{T0_7:.4g}" +
                                             " & " + f"{dmdt_g:.4g}" + " \\\\\n")
                station_properties_table += ("8 & " + f"{p0_8 * 1e-3:.4g}" + " & " + f"{T0_8:.4g}" +
                                             " & " + f"{dmdt_g:.4g}" + " \\\\\n")

                station_properties_table += "\\hline\n"
                station_properties_table += "\\end{tabular}\n"
                station_properties_table += "\\label{tab:thermStatCool_" + f"{coolingFraction:.2g}" + "}\n"
                station_properties_table += "\\end{subtable}"

                overall_performance_table = "\\begin{subtable}{0.45\\textwidth}\n"
                overall_performance_table += "\\centering\n"
                overall_performance_table += "\\caption{Overall Performance}\n"
                overall_performance_table += "\\begin{tabular}{|l|l|}\n"
                overall_performance_table += "\\hline\n"
                overall_performance_table += "Parameter & Value \\\\\n"
                overall_performance_table += "\\hline\n"

                overall_performance_data = [
                    ["Thrust (kN)", f"{Net_thrust * 1e-3:.4g}"],
                    ["Intake mass flow (kg/s)", f"{dmdt_0:.4g}"],
                    ["SFC (mg/Ns)", f"{SFC * 1e6:.4g}"],
                    ["Hot nozzle choke status", "choked" if hotIsChoked else "not choked"],
                    ["Cold nozzle choke status", "choked" if coldIsChoked else "not choked"],
                    ["Hot nozzle pressure ratio", f"{hot_pratio:.4g}"],
                    ["Cold nozzle pressure ratio", f"{cold_pratio:.4g}"],
                    ["Propulsion efficiency", f"{eta_p:.4g}"],
                    ["Thermal efficiency", f"{eta_th:.4g}"],
                    ["Total efficiency", f"{eta_0:.4g}"],
                ]

                for data in overall_performance_data:
                    overall_performance_table += " & ".join(data) + " \\\\\n"

                overall_performance_table += "\\hline\n"
                overall_performance_table += "\\end{tabular}\n"
                overall_performance_table += "\\label{tab:ovPerfCool_" + f"{coolingFraction:.2g}" + "}\n"
                overall_performance_table += "\\end{subtable}\n"
                overall_performance_table += "\\label{tab:anResCool_" + f"{coolingFraction:.2g}" + "}\n"
                overall_performance_table += "\\end{table}\n"

                print(station_properties_table)
                print(overall_performance_table)

        else:
            if hasCooling:
                print('\n\n---Results with cooling flow---')
            else:
                print('\n\n---Results without cooling flow---')

            print('\n--Station thermodynamic properties--')
            print(f'Station 1: P_01: {p0_1 * 1e-3:.4g} kPa, T_01: {T0_1:.4g} K, dmdt_01: {dmdt_0:.3g} kg/s.')
            print(f'Station 2: P_02: {p0_2 * 1e-3:.4g} kPa, T_02: {T0_2:.4g} K, dmdt_02: {dmdt_hot:.3g} kg/s.')
            print(f'Station 3: P_03: {p0_3 * 1e-3:.4g} kPa, T_03: {T0_3:.4g} K, dmdt_03: {dmdt_hot:.3g} kg/s.')
            print(f'Station 4: P_04: {p0_4 * 1e-3:.4g} kPa, T_04: {T0_4:.4g} K, dmdt_04: {dmdt_hot:.3g} kg/s.')
            if hasCooling:
                print(f'Station 5.1: P_05.1: {p0_5 * 1e-3:.4g} kPa, T_05.1: {T0_5:.4g} K, dmdt_05.1: '
                      f'{dmdt_g_1:.3g} kg/s.')
                print(f'Station 5.2: P_05.2: - kPa, T_05.2: {T0_5_2:.4g} K, dmdt_05.2: '
                      f'{dmdt_g_2:.3g} kg/s.')
            else:
                print(f'Station 5: P_05: {p0_5 * 1e-3:.4g} kPa, T_05: {T0_5:.4g} K, dmdt_05: {dmdt_g:.3g} kg/s.')
            print(f'Station 6: P_06: {p0_6 * 1e-3:.4g} kPa, T_06: {T0_6:.4g} K, dmdt_06: {dmdt_g:.3g} kg/s.')
            print(f'Station 7: P_07: {p0_7 * 1e-3:.4g} kPa, T_07: {T0_7:.4g} K, dmdt_07: {dmdt_g:.3g} kg/s.')
            print(f'Station 8: P_08: {p0_8 * 1e-3:.4g} kPa, T_08: {T0_8:.4g} K, dmdt_08: {dmdt_g:.3g} kg/s.')

            print('\n--Overall performance--')
            print(f'Thrust: {Net_thrust * 1e-3:.3g} kN.')
            print(f'Intake mass flow: {dmdt_0:.3g} kg/s.')
            print(f'SFC: {SFC * 1e6:.3g} mg/Ns.')
            print(f'Hot nozzle is {"choked" if hotIsChoked else "not choked"}.')
            print(f'Cold nozzle is {"choked" if coldIsChoked else "not choked"}.')
            print(f'Hot nozzle pressure ratio: {hot_pratio:.3g}')
            print(f'Cold nozzle pressure ratio: {cold_pratio:.3g}')
            print(f'Propulsion efficiency: {eta_p:.3g}.')
            print(f'Thermal efficiency: {eta_th:.3g}.')
            print(f'Total efficiency: {eta_0:.3g}.')

    return F_net, SFC, hotIsChoked, coldIsChoked, hot_pratio, cold_pratio, eta_p, eta_th, eta_0, xi, dmdt, T0, p0, gamma, R, cp


#Assemble mass flow data
massFlowToThrust(1, coolingFraction=0.2, coolSplitFrac=0.4)
massFlows = np.linspace(400, 800, 1000)
testTrust = []
testTrustCool = []

for dmdt in massFlows:
    testTrust.append(massFlowToThrust(dmdt)[0])
    testTrustCool.append(massFlowToThrust(dmdt, coolingFraction=0.2, coolSplitFrac=0.4)[0])
testTrust = np.array(testTrust)
testTrustCool = np.array(testTrustCool)

# Calculate correct mass flow value and print
i_closest = np.argmin(np.abs(testTrust - Net_thrust))
final_dmdt = np.interp(Net_thrust, testTrust[i_closest:i_closest+2], massFlows[i_closest:i_closest+2])

i_closest_c = np.argmin(np.abs(testTrustCool - Net_thrust))
final_dmdt_c = np.interp(Net_thrust, testTrust[i_closest_c:i_closest_c+2], massFlows[i_closest_c:i_closest_c+2])

# Print result from code
pl = False
massFlowToThrust(final_dmdt, hasPrinting=True, printLatex=pl)
dmdt, T0, p0, gamma, R, cp = massFlowToThrust(final_dmdt_c, coolingFraction=0.2, coolSplitFrac=0.4, hasPrinting=True, printLatex=pl)[-6:]

# Plot mass flow to thrust
fig = plt.figure()
plt.hlines(Net_thrust * 1e-3, min(massFlows), max(massFlows),
           label='Thrust requirement', color='g', linestyles='--', zorder=0)
plt.plot(massFlows, testTrust * 1e-3, label='Calculated thrust', color='b', zorder=1)
plt.plot(massFlows, testTrustCool * 1e-3, label='Calculated thrust w. cooling', color='m', zorder=2)
plt.scatter(final_dmdt, Net_thrust * 1e-3, c='r', marker='o', label='Design point', zorder=3)
plt.scatter(final_dmdt_c, Net_thrust * 1e-3, c='k', marker='o', label='Design point w. cooling', zorder=4)
plt.title('Mass flow versus Thrust')
plt.xlabel('Intake mass flow [kg/s]')
plt.ylabel('Net thrust [kN]')
plt.grid()
plt.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/MassFlowToThrust.png', dpi=figureDPI)

# DT1b
#
# BPR_min = 8
# BPR_max = 50
# FPR_min = 1.1
# FPR_max = 1.65
#
# BPR = np.linspace(BPR_min, BPR_max, 100)
# FPR = np.linspace(FPR_min, FPR_max, 100)
#
# sfc = np.zeros((100, 100))
# thrust = np.zeros((100, 100))
# xi = np.zeros((100, 100))
#
# for i, bpr in enumerate(BPR):
#     for j, fpr in enumerate(FPR):
#         thrust[j, i] = massFlowToThrust(1, coolingFraction=0.2, coolSplitFrac=0.4, BPR=bpr, FPR=fpr)[0]
#         massFlow = Net_thrust / thrust[j, i]
#         sfc[j, i] = massFlowToThrust(massFlow, coolingFraction=0.2, coolSplitFrac=0.4, BPR=bpr, FPR=fpr)[1]
#         xi[j, i] = massFlowToThrust(massFlow, coolingFraction=0.2, coolSplitFrac=0.4, BPR=bpr, FPR=fpr)[9]
#         if xi[j, i] > 1.5:
#             xi[j, i] = 1.5
#
# print('hello :D')
#
# plt.figure()
# plt.contourf(BPR, FPR, sfc*1e6, 25)
# plt.colorbar(label='SFC [mg/NS]')
# plt.xlabel('BPR')
# plt.ylabel('FPR')
# # plt.scatter(11.9, 1.55, label='DT1', c='red')
# # plt.vlines(11.9, color='red', ymin=FPR_min, ymax=1.55, linestyle='dashed')
# # plt.hlines(1.55, color='red', xmin=BPR_min, xmax=11.9, linestyle='dashed')
# plt.legend()
#
# plt.figure()
# plt.contourf(BPR, FPR, xi, 25)
# plt.colorbar(label='$\\xi$ [-]')
# plt.xlabel('BPR')
# plt.ylabel('FPR')
# # plt.scatter(11.9, 1.55, label='DT1', c='red')
# # plt.vlines(11.9, color='red', ymin=FPR_min, ymax=1.55, linestyle='dashed')
# # plt.hlines(1.55, color='red', xmin=BPR_min, xmax=11.9, linestyle='dashed')
# plt.legend()
#
# # plt.figure()
# # plt.contourf(BPR, FPR, thrust, 1000)
# # plt.colorbar(label='Thrust [N]')
# # plt.xlabel('BPR')
# # plt.ylabel('FPR')
#
# plt.figure()
# plt.contourf(BPR, FPR, Net_thrust / thrust, 25)
# plt.colorbar(label='mass flow [kg/s]')
# plt.xlabel('BPR')
# plt.ylabel('FPR')
# # plt.scatter(11.9, 1.55, label='DT1', c='red')
# # plt.vlines(11.9, color='red', ymin=FPR_min, ymax=1.55, linestyle='dashed')
# plt.hlines(1.55, color='red', xmin=BPR_min, xmax=11.9, linestyle='dashed')

# DT2 --------------------------------


def areaFuntion(M, station):
    p = p0[station] / (1 + (gamma[station]-1) / 2 * M**2)**gamma[station]/(gamma[station]-1)
    T = T0[station] / (1 + (gamma[station]-1) / 2 * M**2)
    rho = p / R[station] / T

    v = np.sqrt(gamma[station]*R[station]*T) * M
    A = dmdt[station] / rho / v

    return A


def qAreaFunction(M, station):
    return dmdt[station] * np.sqrt(R[station] * T0[station]) / p0[station] / \
        (np.sqrt(gamma[station]) * M * (1 + (gamma[station]-1)/2 * M**2)**(-(gamma[station]+1)/(2*(gamma[station]-1))))


def getRadius(A, htr):
    r_t = np.sqrt(A / np.pi * (1 / (1 - htr ** 2)))

    r_h = htr * r_t

    r_m = (r_t + r_h) / 2

    return r_t, r_h, r_m


def calcStageLoad(dH, U_first, N, r_m_first, r_m_last):
    U_last = r_m_last * U_first / r_m_first
    if N > 1:
        sumSqU = 0
        for i in range(1, N + 1):
            sumSqU += (U_first + (U_last - U_first) * (i - 1) / (N-1)) ** 2
    else:
        sumSqU = ((U_first + U_last) / 2) ** 2

    stageLoad = 2 * dH / sumSqU

    return stageLoad


def calcSOS(M, station):
    T = T0[station] / (1 + (gamma[station]-1) / 2 * M**2)
    return np.sqrt(gamma[station] * R[station] * T)


def calcHubTip(r_mean, area):
    r_tip = area / (4 * np.pi * r_mean) + r_mean
    r_hub = 2 * r_mean - r_tip
    return r_hub, r_tip


def calcAN2(omega, A):
    rps = omega / (2 * np.pi)
    val = A * rps ** 2
    return val


EIS = 2020

# Fan ---------------------------------------------------------------------------
M_ax1_fan = 0.603

htr_1_fan = 44.29 / (98.94 + np.exp(0.0185 * EIS - 33.31))
A_1_fan = qAreaFunction(M_ax1_fan, 1)

r_t1_fan, r_h1_fan = getRadius(A_1_fan, htr_1_fan)[0:2]

r_t2_fan = 0.98 * r_t1_fan

AR_1_fan = 2.4

alpha_fan = 15*np.pi/180    # rad

l_ax1_fan = (1.98*r_t1_fan - 2*r_h1_fan) / (2*AR_1_fan + np.tan(alpha_fan))

r_h2_fan = r_h1_fan + l_ax1_fan * np.tan(alpha_fan)

U_t1_fan = np.sqrt(T0[1]) * (-59.74*FPR + 88.07 * FPR**2 - 25.93 * FPR**3)

omega_fan = U_t1_fan / r_t1_fan

A_2_fan = np.pi*(r_t2_fan**2 - r_h2_fan**2)

A_duct_entry = A_2_fan / (BPR + 1)

r_splitter_lip = np.sqrt(A_duct_entry / np.pi + r_h2_fan**2)

# IPC --------------------------------------------------------------------------------------
AR_duct_FAN_IPC = 0.4

l_ax_duct_FAN_IPC = (r_splitter_lip - r_h2_fan) / AR_duct_FAN_IPC

M_ax_1_IPC = 0.539
M_ax_3_IPC = 0.341

M_t_IPC = 1.55

htr_1_IPC = 0.63
htr_3_IPC = 0.819

psi_IPC = -8.968 + 0.004877 * EIS

A_1_IPC = qAreaFunction(M_ax_1_IPC, 2)
A_3_IPC = qAreaFunction(M_ax_3_IPC, 3)

r_t1_IPC, r_h1_IPC, r_m1_IPC = getRadius(A_1_IPC, htr_1_IPC)
r_t3_IPC, r_h3_IPC, r_m3_IPC = getRadius(A_3_IPC, htr_3_IPC)

dH_IPC = cp[2] * (T0[3] - T0[2])

# Interpolate to get requested blade tip speed
C_t1_IPC_rel = M_t_IPC * calcSOS(M_ax_1_IPC, 2)
C_ax1_IPC = M_ax_1_IPC * calcSOS(M_ax_1_IPC, 2)  # Assume no swirl
U_t1_IPC = np.sqrt(C_t1_IPC_rel ** 2 - C_ax1_IPC ** 2)

omega_IPC = U_t1_IPC / r_t1_IPC

U_m1_IPC = r_m1_IPC * omega_IPC
U_m3_IPC = r_m3_IPC * omega_IPC
U_t3_IPC = r_t3_IPC * omega_IPC
U_h1_IPC = r_h1_IPC * omega_IPC
U_h3_IPC = r_h3_IPC * omega_IPC

idealres = []
N_list = []
for i in range(2, 8):
    N_list.append(i)
    idealres.append(np.abs(calcStageLoad(dH_IPC, U_m1_IPC, i, r_m1_IPC, r_m3_IPC) - psi_IPC))

N_stages_IPC = N_list[np.argmin(idealres)]

# aspect ratios IPC
AR_1_IPC = 36.20 - 0.01694*EIS
AR_3_IPC = 35.47 - 0.01694*EIS

h_1_IPC = r_t1_IPC - r_h1_IPC
h_3_IPC = r_t3_IPC - r_h3_IPC
h_mean_1_IPC = np.sqrt(h_1_IPC*h_3_IPC)
AR_mean_IPC = np.sqrt(AR_1_IPC * AR_3_IPC)
c = 0.3  # spacing
l_ax_IPC = 2 * N_stages_IPC * h_mean_1_IPC * (1 + c) / AR_mean_IPC

# HPC --------------------------------------------------------------------
AR_duct_IPC_HPC = 0.4

M_ax_1_HPC = 0.482
M_ax_3_HPC = 0.263

M_t_HPC = 1.3

htr_1_HPC = 0.5613 / (1.487 - np.exp(-0.04286 * (p0[4]/p0[3] + 0.5718)))
htr_3_HPC = 0.908

psi_HPC = -5.736 + 0.00323 * EIS

A_1_HPC = qAreaFunction(M_ax_1_HPC, 3)
A_3_HPC = qAreaFunction(M_ax_3_HPC, 4)

r_t1_HPC, r_h1_HPC, r_m1_HPC = getRadius(A_1_HPC, htr_1_HPC)
r_t3_HPC, r_h3_HPC, r_m3_HPC = getRadius(A_3_HPC, htr_3_HPC)

dH_HPC = cp[3] * (T0[4] - T0[3])

# Interpolate to get requested blade tip speed
C_t1_HPC_rel = M_t_HPC * calcSOS(M_ax_1_HPC, 3)
C_ax1_HPC = M_ax_1_HPC * calcSOS(M_ax_1_HPC, 3)  # Assume no swirl
U_t1_HPC = np.sqrt(C_t1_HPC_rel ** 2 - C_ax1_HPC ** 2)

omega_HPC = U_t1_HPC / r_t1_HPC

U_m1_HPC = r_m1_HPC * omega_HPC
U_m3_HPC = r_m3_HPC * omega_HPC
U_t3_HPC = r_t3_HPC * omega_HPC
U_h1_HPC = r_h1_HPC * omega_HPC
U_h3_HPC = r_h3_HPC * omega_HPC

idealres = []
N_list = []
for i in range(2, 8):
    N_list.append(i)
    idealres.append(np.abs(calcStageLoad(dH_HPC, U_m1_HPC, i, r_m1_HPC, r_m3_HPC) - psi_HPC))

N_stages_HPC = N_list[np.argmin(idealres)]
psi_HPC = calcStageLoad(dH_HPC, U_m1_HPC, N_stages_HPC, r_m1_HPC, r_m3_HPC)

# aspect ratios IPC
AR_1_HPC = 31.40 - 0.0147*EIS
AR_3_HPC = 30.70 - 0.0147*EIS

h_1_HPC = r_t1_HPC - r_h1_HPC
h_3_HPC = r_t3_HPC - r_h3_HPC
h_mean_1_HPC = np.sqrt(h_1_HPC*h_3_HPC)
AR_mean_HPC = np.sqrt(AR_1_HPC * AR_3_HPC)
c = 0.3  # spacing
l_ax_HPC = 2 * N_stages_HPC * h_mean_1_HPC * (1 + c) / AR_mean_HPC

l_ax_duct_IPC_HPC = (r_t3_IPC - r_h3_IPC + r_t1_HPC - r_h1_HPC) / (2 * AR_duct_IPC_HPC)

# Combustor -------------------------------------------------------------------------------------
t_CC = 6e-3  # s combustor residence time
M_CC = 0.06  # - Mach number average in combustor
V_CC = M_CC * calcSOS(M_CC, 4)

l_CC = t_CC * V_CC

beta_CC = 10  # degrees
dy_CC = l_CC * np.tan(np.pi / 180 * beta_CC)

r_m1_HPT = r_m3_HPC + dy_CC

# HPT ------------------------------------------------------------------------------------------
psi_HPT_required = 3.247
omega_HPT = omega_HPC

dH_HPT = cp[5] * (T0[5] + (1562 - 1680) - T0[6])
U_m1_HPT = r_m1_HPT * omega_HPC

M_ax_1_HPT = 0.150
M_ax_3_HPT = 0.420

A_1_HPT = qAreaFunction(M_ax_1_HPT, 5)
A_3_HPT = qAreaFunction(M_ax_3_HPT, 6)

r_m3_HPT = r_m1_HPT  # Assume same radius

r_h1_HPT, r_t1_HPT = calcHubTip(r_m1_HPT, A_1_HPT)
r_h3_HPT, r_t3_HPT = calcHubTip(r_m3_HPT, A_3_HPT)

# Find if number of stages is one or two
N_stages_HPT = 1
psi_HPT_single_stage = 2 * dH_HPT / (U_m1_HPT ** 2)
psi_HPT_double_stage = 2 * dH_HPT / (2 * U_m1_HPT ** 2)

if psi_HPT_single_stage > psi_HPT_required:
    N_stages_HPT = 2

if psi_HPT_double_stage > psi_HPT_required:
    print('Error: Number of hpt stages exceeds 2!')
    exit()

AN2_HPT = calcAN2(omega_HPT, A_1_HPT)

if N_stages_HPT == 1:
    AN2crit_HPT = (-158.99 + EIS * 0.0830) * 1e3
    if AN2_HPT > AN2crit_HPT:
        print('Error: Exceeding stress predictions in HPT!')
        exit()

elif N_stages_HPT == 2:
    AN2crit_1_HPT = (-129.44 + EIS * 0.0675) * 1e3
    AN2crit_3_HPT = (-177.33 + EIS * 0.0930) * 1e3

    if AN2_HPT > AN2crit_1_HPT or AN2_HPT > AN2crit_3_HPT:
        print('Error: Exceeding stress predictions in HPT!')
        exit()

# Calculate axial length
if N_stages_HPT == 1:
    AR_1_HPT = 29.233 - 0.0140 * EIS

    h_1_HPT = r_t1_HPT - r_h1_HPT
    h_3_HPT = r_t3_HPT - r_h3_HPT
    h_mean_1_HPT = np.sqrt(h_1_HPT * h_3_HPT)
    c = 0.2  # spacing
    l_ax_HPT = 2 * N_stages_HPT * h_mean_1_HPT * (1 + c) / AR_1_HPT

else:
    AR_1_HPT = 29.510 - 0.0140 * EIS
    AR_3_HPT = 30.143 - 0.0140 * EIS

    h_1_HPT = r_t1_HPT - r_h1_HPT
    h_3_HPT = r_t3_HPT - r_h3_HPT
    h_mean_1_HPT = np.sqrt(h_1_HPT * h_3_HPT)
    AR_mean_HPT = np.sqrt(AR_1_HPT * AR_3_HPT)
    c = 0.2  # spacing
    l_ax_HPT = 2 * N_stages_HPT * h_mean_1_HPT * (1 + c) / AR_mean_HPT


# IPT ------------------------------------------------------------------------------------------
psi_IPT_required = 3.247
omega_IPT = omega_IPC

r_m1_IPT = 1.3 * r_m3_HPT
r_m3_IPT = r_m1_IPT

dH_IPT = cp[6] * (T0[6] - T0[7])
U_m1_IPT = r_m1_IPT * omega_IPT

M_ax_1_IPT = 0.314
M_ax_3_IPT = 0.554

A_1_IPT = qAreaFunction(M_ax_1_IPT, 6)
A_3_IPT = qAreaFunction(M_ax_3_IPT, 7)

r_h1_IPT, r_t1_IPT = calcHubTip(r_m1_IPT, A_1_IPT)
r_h3_IPT, r_t3_IPT = calcHubTip(r_m3_IPT, A_3_IPT)

# Find if number of stages is one or two
N_stages_IPT = 1
psi_IPT_single_stage = 2 * dH_IPT / (U_m1_IPT ** 2)
psi_IPT_double_stage = 2 * dH_IPT / (2 * U_m1_IPT ** 2)

if psi_IPT_single_stage > psi_IPT_required:
    N_stages_IPT = 2

if psi_IPT_double_stage > psi_IPT_required:
    print('Error: Number of hpt stages exceeds 2!')
    exit()

AN2_IPT = calcAN2(omega_IPT, A_1_IPT)
AN2crit_1_IPT = (-133.97 + EIS * 0.07099) * 1e3

if AN2_IPT > AN2crit_1_IPT:
    print('Error: Exceeding stress predictions in IPT!')
    exit()

# Calculate axial length
AR_1_IPT = 30.492 - 0.0140 * EIS

h_1_IPT = r_t1_IPT - r_h1_IPT
h_3_IPT = r_t3_IPT - r_h3_IPT
h_mean_1_IPT = np.sqrt(h_1_IPT * h_3_IPT)
c = 0.2  # spacing
l_ax_IPT = 2 * N_stages_IPT * h_mean_1_IPT * (1 + c) / AR_1_IPT

# LPT -------------------------------------------------------------------------------------
r_m1_LPT = 1.15 * r_m3_IPT
omega_LPT = omega_fan
U_m1_LPT = r_m1_LPT * omega_LPT

ttr_LPT_FAN = 1.25 / (1.9176 - np.exp(-0.2503 * (BPR + 0.6410)))
r_t3_LPT = r_t1_fan * ttr_LPT_FAN

psi_LPT = -39.26 + 0.02185 * EIS

M_ax_1_LPT = 0.368
M_ax_3_LPT = 0.322

A_1_LPT = qAreaFunction(M_ax_1_LPT, 7)
A_3_LPT = qAreaFunction(M_ax_3_LPT, 8)

r_h1_LPT, r_t1_LPT = calcHubTip(r_m1_LPT, A_1_LPT)
r_h3_LPT = np.sqrt(r_t3_LPT ** 2 - A_3_LPT / np.pi)
r_m3_LPT = (r_t3_LPT + r_h3_LPT) / 2

AR_1_LPT = 31.40 - 0.0146 * EIS
AR_3_LPT = -35.60 + 0.0209 * EIS

dH_LPT = cp[7] * (T0[7] - T0[8])

idealres = []
N_list = []
for i in range(2, 20):
    N_list.append(i)
    idealres.append(np.abs(calcStageLoad(dH_LPT, U_m1_LPT, i, r_m1_LPT, r_m3_LPT) - psi_LPT))

N_stages_LPT = N_list[np.argmin(idealres)]

h_1_LPT = r_t1_LPT - r_h1_LPT
h_3_LPT = r_t3_LPT - r_h3_LPT
h_mean_1_LPT = np.sqrt(h_1_LPT * h_3_LPT)
AR_mean_LPT = np.sqrt(AR_1_LPT * AR_3_LPT)
c = 0.4  # spacing
l_ax_LPT = 2 * N_stages_LPT * h_mean_1_LPT * (1 + c) / AR_mean_LPT

AN2_1_LPT = calcAN2(omega_LPT, A_1_LPT)
AN2_3_LPT = calcAN2(omega_LPT, A_3_LPT)

# DUCT IPT-LPT--------------------------
h_mean_duct_HPT_IPT = (r_t3_HPT - r_h3_HPT) + (r_t1_IPT - r_h1_IPT)
l_ax_duct_HPT_IPT = h_mean_duct_HPT_IPT / 0.4

# DUCT IPT-LPT--------------------------
h_mean_duct_IPT_LPT = (r_t3_IPT - r_h3_IPT) + (r_t1_LPT - r_h1_LPT)
l_ax_duct_IPT_LPT = h_mean_duct_IPT_LPT / 0.4

# drawing engine sketch ------------------------------------------------------------------------
# ul = upper left
# ur = upper right
# ll = lower left
# lr = lower right
# [ul, ur, lr, ll, ul]

# fan
x_fan = [0, l_ax1_fan, l_ax1_fan, 0, 0]
y_fan = [r_t1_fan, r_t2_fan, r_h2_fan, r_h1_fan, r_t1_fan]

# IPC
x0_IPC = l_ax1_fan + l_ax_duct_FAN_IPC

x_IPC = [x0_IPC, x0_IPC + l_ax_IPC, x0_IPC + l_ax_IPC, x0_IPC, x0_IPC]
y_IPC = [r_t1_IPC, r_t3_IPC, r_h3_IPC, r_h1_IPC, r_t1_IPC]

# HPC
x0_HPC = x0_IPC + + l_ax_IPC + l_ax_duct_IPC_HPC

x_HPC = [x0_HPC, x0_HPC + l_ax_HPC, x0_HPC + l_ax_HPC, x0_HPC, x0_HPC]
y_HPC = [r_t1_HPC, r_t3_HPC, r_h3_HPC, r_h1_HPC, r_t1_HPC]

# CC
x0_CC = x0_HPC + l_ax_HPC

x_CC = [x0_CC, x0_CC + l_CC, x0_CC + l_CC, x0_CC, x0_CC]
y_CC = [r_t3_HPC, r_t1_HPT, r_h1_HPT, r_h3_HPC, r_t3_HPC]

# HPT
x0_HPT = x0_CC + l_CC

x_HPT = [x0_HPT, x0_HPT + l_ax_HPT, x0_HPT + l_ax_HPT, x0_HPT, x0_HPT]
y_HPT = [r_t1_HPT, r_t3_HPT, r_h3_HPT, r_h1_HPT, r_t1_HPT]

# IPT
x0_IPT = x0_HPT + l_ax_HPT + l_ax_duct_HPT_IPT

x_IPT = [x0_IPT, x0_IPT + l_ax_IPT, x0_IPT + l_ax_IPT, x0_IPT, x0_IPT]
y_IPT = [r_t1_IPT, r_t3_IPT, r_h3_IPT, r_h1_IPT, r_t1_IPT]

#LPT
x0_LPT = x0_IPT + l_ax_IPT + l_ax_duct_IPT_LPT

x_LPT = [x0_LPT, x0_LPT + l_ax_LPT, x0_LPT + l_ax_LPT, x0_LPT, x0_LPT]
y_LPT = [r_t1_LPT, r_t3_LPT, r_h3_LPT, r_h1_LPT, r_t1_LPT]

# Plot engine layout
fig, ax = plt.subplots()
ax.axhline(y=0, color='k', linestyle='--')
ax.fill(x_fan, y_fan, color='b', label='Fan/low pressure')
ax.fill(x_IPC, y_IPC, color='g', label='Intermediate pressure')
ax.fill(x_HPC, y_HPC, color='c', label='High pressure')
ax.fill(x_CC, y_CC, color='r', label='Compustion chamber')
ax.fill(x_HPT, y_HPT, color='c')
ax.fill(x_IPT, y_IPT, color='g')
ax.fill(x_LPT, y_LPT, color='b')
ax.set_ylim(-0.1, r_t1_fan+0.1)
ax.axis('equal')
ax.set_xlabel('x [m]')
ax.set_ylabel('r [m]')
ax.set_title('Engine component layout')
ax.grid()
ax.legend()

# Save figure
figureDPI = 200
fig.set_size_inches(8, 6)
fig.savefig('img/EngineLayout.png', dpi=figureDPI)

# prints for table
print('\n--Engine sizing---')
printLatex = True
if printLatex:
    component_lengths_table = "\\begin{table}[ht]\n"
    component_lengths_table += "\\centering\n"
    component_lengths_table += "\\caption{Engine sizing results}\n"
    component_lengths_table += "\\begin{subtable}{0.45\\textwidth}\n"
    component_lengths_table += "\\centering\n"
    component_lengths_table += "\\caption{Component Lengths}\n"
    component_lengths_table += "\\begin{tabular}{|l|l|}\n"
    component_lengths_table += "\\hline\n"
    component_lengths_table += "Component & Length (m) \\\\\n"
    component_lengths_table += "\\hline\n"

    component_lengths_data = [
        ["Fan", f"{l_ax1_fan:.4g}"],
        ["IPC", f"{l_ax_IPC:.4g}"],
        ["HPC", f"{l_ax_HPC:.4g}"],
        ["Combustion Chamber", f"{l_CC:.4g}"],
        ["HPT", f"{l_ax_HPT:.4g}"],
        ["IPT", f"{l_ax_IPT:.4g}"],
        ["LPT", f"{l_ax_LPT:.4g}"],
        ["Total Length", f"{x0_LPT + l_ax_LPT:.4g}"],
    ]

    for data in component_lengths_data:
        component_lengths_table += " & ".join(data) + " \\\\\n"

    component_lengths_table += "\\hline\n"
    component_lengths_table += "\\end{tabular}\n"
    component_lengths_table += "\\label{tab:compLength}\n"
    component_lengths_table += "\\end{subtable}\n"

    num_of_stages_table = "\\begin{subtable}{0.45\\textwidth}\n"
    num_of_stages_table += "\\centering\n"
    num_of_stages_table += "\\caption{Number of Stages}\n"
    num_of_stages_table += "\\begin{tabular}{|l|l|}\n"
    num_of_stages_table += "\\hline\n"
    num_of_stages_table += "Component & Number of Stages \\\\\n"
    num_of_stages_table += "\\hline\n"

    num_of_stages_data = [
        ["Fan", "1"],
        ["IPC", str(N_stages_IPC)],
        ["HPC", str(N_stages_HPC)],
        ["HPT", str(N_stages_HPT)],
        ["IPT", str(N_stages_IPT)],
        ["LPT", str(N_stages_LPT)],
    ]

    for data in num_of_stages_data:
        num_of_stages_table += " & ".join(data) + " \\\\\n"

    num_of_stages_table += "\\hline\n"
    num_of_stages_table += "\\end{tabular}\n"
    num_of_stages_table += "\\label{tab:stageNo}\n"
    num_of_stages_table += "\\end{subtable}\n"
    num_of_stages_table += "\\label{tab:engSizing}\n"
    num_of_stages_table += "\\end{table}"

    print(component_lengths_table)
    print(num_of_stages_table)
else:
    print('Component lengths:')
    print(f'Fan: {l_ax1_fan: .4g} m')
    print(f'IPC: {l_ax_IPC: .4g} m')
    print(f'HPC: {l_ax_HPC: .4g} m')
    print(f'Combustion chamber: {l_CC: .4g} m')
    print(f'HPT: {l_ax_HPT: .4g} m')
    print(f'IPT: {l_ax_IPT: .4g} m')
    print(f'LPT: {l_ax_LPT: .4g} m')
    print(f'Total length: {x0_LPT + l_ax_LPT: .4g} m')

    print('Number of stages:')
    print(f'Fan 1')
    print(f'IPC {N_stages_IPC}')
    print(f'HPC {N_stages_HPC}')
    print(f'HPT {N_stages_HPT}')
    print(f'IPT {N_stages_IPT}')
    print(f'LPT {N_stages_LPT}')

print('\n--Turbine stresses (AN2)--')
if printLatex:
    single_table_package = "\\begin{table}[ht]\n"
    single_table_package += "\\centering\n"
    single_table_package += "\\begin{tabular}{|l|l|}\n"
    single_table_package += "\\hline\n"
    single_table_package += "Component & Stress (m\\textsuperscript{2}/s\\textsuperscript{2}) \\\\\n"
    single_table_package += "\\hline\n"

    data = [
        ["$AN^2_\\text{HPT}$", f"{AN2_HPT:.4g}"],
        ["$AN^2_\\text{IPT}$", f"{AN2_IPT:.4g}"],
        ["$AN^2_\\text{LPT, 1}$", f"{AN2_1_LPT:.4g}"],
        ["$AN^2_\\text{LPT, 3}$", f"{AN2_3_LPT:.4g}"],
    ]

    for entry in data:
        single_table_package += " & ".join(entry) + " \\\\\n"

    single_table_package += "\\hline\n"
    single_table_package += "\\end{tabular}\n"
    single_table_package += "\\caption{AN2 Parameters}\n"
    single_table_package += "\\end{table}\n"

    print(single_table_package)
else:
    print(f'AN2_HPT = {AN2_HPT:.3g} m^2 / s^2')
    print(f'AN2_IPT = {AN2_IPT:.3g} m^2 / s^2')
    print(f'AN2_1_LPT = {AN2_1_LPT:.3g} m^2 / s^2')
    print(f'AN2_3_LPT = {AN2_3_LPT:.3g} m^2 / s^2')

# DT3 --------------------------------------------

# Calculate velocity and sos at HPC entry
a_1_HPC = calcSOS(M_ax_1_HPC, 3)
v_1_HPC = a_1_HPC * M_ax_1_HPC  # Assume no swirl here

# Calculate velocity and SOS at HPC first stage exit
lenFracHPC = (2 + 0.3) * h_1_HPC / AR_1_HPC / l_ax_HPC
M_ax_2_HPC = M_ax_1_HPC * (1 + lenFracHPC * (M_ax_3_HPC / M_ax_1_HPC - 1))

a = 11
plt.show()
