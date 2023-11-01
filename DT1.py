import numpy as np


# Given parameters
Altitude = 35000 * 0.3048  # m
Flight_Mach_number = 0.82
dT_isa = 10  # K
Net_thrust = 67.3 * 1e3  # N
BPR = 11.8
FPR = 1.55
OPR = 47
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

# Calculate speed of sound at design altitude
sourceAlt = [10000, 11000]
sourceSOS = [295.2, 295.1]
Speed_of_sound = np.interp(Altitude, sourceAlt, sourceSOS)  # m/s

C_a = Flight_Mach_number * Speed_of_sound
a = 10
