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


def massFlowToThrust(dmdt_0, coolingFraction=0.0, coolSplitFrac=0.0, hasPrinting=False, printLatex=False):
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
    dWdt_fan = dmdt_hot * cp_a * (T0_2 - T0_1)

    p0_3 = p0_2 * IPC_pressure_ratio    # pressure after IPC
    T0_3 = T0_2 * IPC_pressure_ratio ** ((gamma_a - 1) / (IPC_polytropic_efficiency * gamma_a))
    dWdt_IPC = dmdt_hot * cp_a * (T0_3 - T0_2)

    p0_4 = p0_3 * HPC_pressure_ratio    # pressure after HPC
    T0_4 = T0_3 * HPC_pressure_ratio ** ((gamma_a - 1) / (HPC_polytropic_efficiency * gamma_a))
    dWdt_HPC = dmdt_hot * cp_a * (T0_4 - T0_3)

    # carlos code
    # TODO static conditions not total
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
    T0_5_2 = (dmdt_CC * cp_g * T0_5 + dmdt_cool_s * cp_a * T0_4) / (cp_g * dmdt_g_2)

    T0_6 = T0_5_2 - dWdt_HPC / (dmdt_g_2 * cp_g * Shaft_mechanical_efficiency)
    p0_6 = p0_5 * (T0_6 / T0_5_2) ** (gamma_g / ((gamma_g - 1) * HPT_polytropic_efficiency))
    dmdt_g = dmdt_g_2 + dmdt_cool_r

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
        T_9 = T0_8-Hot_jet_efficiency*T0_8*(1-(1/p0_8/p_9)**((gamma_g-1)/gamma_g))
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
        cold_pratio = CPR_hot
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

    if hasPrinting:
        if printLatex:
            if hasPrinting:
                station_properties_table = "\n\\begin{table}[ht]\n"
                station_properties_table += "\\centering\n"
                station_properties_table += "\\begin{tabular}{|l|l|l|l|}\n"
                station_properties_table += "\\hline\n"
                station_properties_table += "Station & Pressure (kPa) & Temperature (K) & Mass Flow (kg/s) \\\\\n"
                station_properties_table += "\\hline\n"

                station_data = [
                    ["1", f"{p0_1 * 1e-3:.4g}", f"{T0_1:.4g}", f"{dmdt_0:.3g}"],
                    ["2", f"{p0_2 * 1e-3:.4g}", f"{T0_2:.4g}", f"{dmdt_hot:.3g}"],
                    ["3", f"{p0_3 * 1e-3:.4g}", f"{T0_3:.4g}", f"{dmdt_hot:.3g}"],
                    ["4", f"{p0_4 * 1e-3:.4g}", f"{T0_4:.4g}", f"{dmdt_hot:.3g}"],
                ]

                for data in station_data:
                    station_properties_table += " & ".join(data) + " \\\\\n"

                if hasCooling:
                    station_data_cooling = [
                        ["5.1", f"{p0_5 * 1e-3:.4g}", f"{T0_5:.4g}", f"{dmdt_g_1:.3g}"],
                        ["5.2", f"{p0_5 * 1e-3:.4g}", f"{T0_5_2:.4g}", f"{dmdt_g_2:.3g}"],
                    ]

                    for data in station_data_cooling:
                        station_properties_table += " & ".join(data) + " \\\\\n"
                else:
                    station_data_no_cooling = [
                        ["5", f"{p0_5 * 1e-3:.4g}", f"{T0_5:.4g}", f"{dmdt_g:.3g}"],
                    ]

                    for data in station_data_no_cooling:
                        station_properties_table += " & ".join(data) + " \\\\\n"

                station_properties_table += "6 & " + f"{p0_6 * 1e-3:.4g}" + " & " + f"{T0_6:.4g}" + " & " + f"{dmdt_g:.3g}" + " \\\\\n"
                station_properties_table += "7 & " + f"{p0_7 * 1e-3:.4g}" + " & " + f"{T0_7:.4g}" + " & " + f"{dmdt_g:.3g}" + " \\\\\n"
                station_properties_table += "8 & " + f"{p0_8 * 1e-3:.4g}" + " & " + f"{T0_8:.4g}" + " & " + f"{dmdt_g:.3g}" + " \\\\\n"

                station_properties_table += "\\hline\n"
                station_properties_table += "\\end{tabular}\n"
                station_properties_table += "\\caption{Station Thermodynamic Properties}\n"
                station_properties_table += "\\end{table}"

                print(station_properties_table)

                overall_performance_table = "\n\\begin{table}[ht]\n"
                overall_performance_table += "\\centering\n"
                overall_performance_table += "\\begin{tabular}{|l|l|}\n"
                overall_performance_table += "\\hline\n"
                overall_performance_table += "Parameter & Value \\\\\n"
                overall_performance_table += "\\hline\n"

                overall_performance_data = [
                    ["Thrust (kN)", f"{Net_thrust * 1e-3:.3g}"],
                    ["Intake mass flow (kg/s)", f"{dmdt_0:.3g}"],
                    ["SFC (mg/Ns)", f"{SFC * 1e6:.3g}"],
                    ["Hot channel choke status", "choked" if hotIsChoked else "not choked"],
                    ["Cold channel choke status", "choked" if coldIsChoked else "not choked"],
                    ["Hot nozzle pressure ratio", f"{hot_pratio:.3g}"],
                    ["Cold nozzle pressure ratio", f"{cold_pratio:.3g}"],
                    ["Propulsion efficiency", f"{eta_p:.3g}"],
                    ["Thermal efficiency", f"{eta_th:.3g}"],
                    ["Total efficiency", f"{eta_0:.3g}"],
                ]

                for data in overall_performance_data:
                    overall_performance_table += " & ".join(data) + " \\\\\n"

                overall_performance_table += "\\hline\n"
                overall_performance_table += "\\end{tabular}\n"
                overall_performance_table += "\\caption{Overall Performance}\n"
                overall_performance_table += "\\end{table}"

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
                print(f'Station 5.2: P_05.2: {p0_5 * 1e-3:.4g} kPa, T_05.2: {T0_5_2:.4g} K, dmdt_05.2: '
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
            print(f'Hot channel is {"choked" if hotIsChoked else "not choked"}.')
            print(f'Cold channel is {"choked" if coldIsChoked else "not choked"}.')
            print(f'Hot nozzle pressure ratio: {hot_pratio:.3g}')
            print(f'Cold nozzle pressure ratio: {cold_pratio:.3g}')
            print(f'Propulsion efficiency: {eta_p:.3g}.')
            print(f'Thermal efficiency: {eta_th:.3g}.')
            print(f'Total efficiency: {eta_0:.3g}.')

    return F_net, SFC, hotIsChoked, coldIsChoked, hot_pratio, cold_pratio, eta_p, eta_th, eta_0


# Assemble mass flow data
massFlows = np.linspace(300, 500, 1000)
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
pl = True
massFlowToThrust(final_dmdt, hasPrinting=True, printLatex=pl)
massFlowToThrust(final_dmdt_c, coolingFraction=0.2, coolSplitFrac=0.4, hasPrinting=True, printLatex=pl)

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

plt.show()
