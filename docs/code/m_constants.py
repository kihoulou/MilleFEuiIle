import numpy as np


# Time related constants
yr  = 3600.0*24*365.25
kyr = 3600.0*24*365.25*1e3
Myr = 3600.0*24*365.25*1e6

# Molar gass constant
R_gas   = 8.314   # J K^-1 mol^-1

# Insolation at 1 AU from sun
insol_1AU = 1360.0*np.sin(np.pi*20.0/180)  # W m^2

# Stefan-Boltzmann constant
SB_const = 5.67e-8 		# W m^−2 K^−4.

# Ice activation volume
V_act = -1.3e-5

