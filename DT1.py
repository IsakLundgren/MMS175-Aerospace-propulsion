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
T_a = np.interp(Altitude, sourceAlt, sourceTemperature)  # K
p_a = np.interp(Altitude, sourceAlt, sourcePressure)  # Pa
Speed_of_sound = np.interp(Altitude, sourceAlt, sourceSOS)  # m/s
gamma = 1.4

# Flight velocity
C_a = Flight_Mach_number * Speed_of_sound  # m/s

# IPC pressure ratio
IPC_pressure_ratio = OPR / (FPR * HPC_pressure_ratio)

# Follow thermodynamic quantities along flow path
# Note that the measurements take place directly after the component
# Assume calorically perfect gas
# TODO Add losses from efficiencies
T0_a = T_a * (1 + (gamma - 1) / 2 * Flight_Mach_number ** 2)
p0_a = p_a * (1 + (gamma - 1) / 2 * Flight_Mach_number ** 2) ** (gamma / (gamma - 1))

p0_fan = p0_a * FPR

p0_IPC = p0_fan * IPC_pressure_ratio

p0_HPC = p0_IPC * HPC_pressure_ratio

# Intake mass flow TODO Complete calculations of prerequisites
dmdt_hot = 1  # kg/s
dmdt_cold = dmdt_hot * BPR  # kg/s
C_hot = 1  # m/s
C_cold = 1  # m/s
A_hot = 1  # m2
A_cold = 1  # m2
p_hot = 1  # Pa
p_cold = 1  # Pa

dmdt_0 = 1 / C_a * (dmdt_hot * C_hot + (p_hot - p_a) * A_hot +
                    dmdt_cold * C_cold + (p_cold - p_a) * A_cold -
                    Net_thrust)

print(f'Intake mass flow: {dmdt_0:.3g} kg/s.')
