import numpy as np


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
Intake_efficiency = 99.5  # %
Fan_polytropic_efficiency = 92.0  # %
Combustor_efficiency = 99.9  # %
Combustor_pressure_loss = 3.5  # %
HPT_polytropic_efficiency = 90.5  # %
IPT_polytropic_efficiency = 91.0  # %
LPT_polytropic_efficiency = 91.5  # %
Cold_Jet_efficiency = 98.0  # %
Hot_jet_efficiency = 99.0  # %
Shaft_mechanical_efficiency = 99.5  # %
IPC_polytropic_efficiency = 91.0  # %
HPC_polytropic_efficiency = 91.5  # %

# Thermodynamic ambient properties at design altitude
sourceAlt = [10000, 11000]  # m https://www.engineeringtoolbox.com/elevation-speed-sound-air-d_1534.html
sourceTemperature = [-49.9 + 273.15, -56.4 + 273.15]  # K
sourcePressure = [26.48 * 1e3, 22.68 * 1e3]  # Pa
sourceSOS = [299.5, 295.2]  # m/s
T_1 = np.interp(Altitude, sourceAlt, sourceTemperature) + dT_isa   # K
p_1 = np.interp(Altitude, sourceAlt, sourcePressure)  # Pa
gamma_a = 1.4
gamma_g = 1.333
cp_a = 1005 # J/kgK
cp_g = 1148 # # J/kgK
R_a = 287     # unit
Speed_of_sound_1 = np.sqrt(gamma_a * R_a * T_1)  # m/s


# Flight velocity
C_a = Flight_Mach_number * Speed_of_sound_1  # m/s

# IPC pressure ratio
IPC_pressure_ratio = OPR / (FPR * HPC_pressure_ratio)


def massFlowToThrust(dmdt_0):
    dmdt_hot = dmdt_0 / (BPR + 1)

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

    dmdt_f = FAR * dmdt_hot
    dmdt_g = dmdt_f + dmdt_hot

    T0_5 = Turbine_inlet_temperature
    p0_5 = p0_4 * Combustor_pressure_loss

    T0_6 = T0_5 - dWdt_HPC / (dmdt_g * cp_g * Shaft_mechanical_efficiency)
    p0_6 = p0_5 * (T0_6 / T0_5) ** (gamma_g / ((gamma_g - 1) * HPT_polytropic_efficiency))

    T0_7 = T0_6 - dWdt_IPC / (dmdt_g * cp_g * Shaft_mechanical_efficiency)
    p0_7 = p0_6 * (T0_7 / T0_6) ** (gamma_g / ((gamma_g - 1) * IPT_polytropic_efficiency))

    T0_8 = T0_7 - dWdt_fan / (dmdt_g * cp_g * Shaft_mechanical_efficiency)
    p0_8 = p0_7 * (T0_8 / T0_7) ** (gamma_g / ((gamma_g - 1) * LPT_polytropic_efficiency))

    # pressure ratio critical
    CPR_hot = (1 / (1 - 1 / Hot_jet_efficiency * (gamma_g - 1) / (gamma_g + 1)) ** (gamma_g/(gamma_g-1)))
    CPR_cold = (1 / (1 - 1 / Cold_Jet_efficiency * (gamma_a - 1) / (gamma_a + 1)) ** (gamma_a/(gamma_a-1)))

    hotIsChoked = (p0_8 / p_1) > po9_pc
    coldIsChoked = (p0_2 / p_1) > p02_p10c

    # TODO if statement for choked and not-choked conditions

    if hotIsChoked: # p9=p9c critical pressure
        p_9= p0_8 / CPR_hot
    else:
        p_9=p_1

    if coldIsChoked:
        p_10 = p0_2 / CPR_cold
    else:
        p_10 = p_1

    # Intake mass flow TODO Complete calculations of prerequisites
    dmdt_hot = 1  # kg/s
    dmdt_cold = dmdt_hot * BPR  # kg/s
    C_hot = 1  # m/s
    C_cold = 1  # m/s
    A_hot = 1  # m2
    A_cold = 1  # m2
    p_hot = 1  # Pa
    p_cold = 1  # Pa

    return (dmdt_hot * C_hot + (p_hot - p0_1) * A_hot +
                    dmdt_cold * C_cold + (p_cold - p0_1) * A_cold -
                    dmdt_0 * C_a), dmdt_0

# Guess loop
first_dmdt_0_guess = 10  # kg/s
# TODO make loop
final_thrust, final_dmdt = massFlowToThrust(first_dmdt_0_guess)

print(f'Intake mass flow: {final_dmdt:.3g} kg/s.')
print(f'Thrust: {final_thrust:.3g} N.')


