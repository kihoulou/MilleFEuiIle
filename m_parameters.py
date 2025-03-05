from dolfin import *
from m_constants import *
import numpy as np

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

# --- Structure of this file ---

# --- RUNNING PARAMETERS ---
# 1/ Output files settings
# 2/ Reloading options

# --- GEOMETRY ---
# 3/ Mesh geometry
# 4/ Material composition

# --- PHYSICS ---
# 5/ Stokes problem settings
# 6/ Heat transport and melting settings
# 7/ Rheology

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 1/ OUTPUT FILES SETTINGS -------------------
#----------------------------------------------------------------------

# --- Name of the directory with results ---
name = "van_Keken"

# --- In what units the time will be? ---
# 1.0 - seconds or nondimensional, or yr, kyr or Myr
time_units = 1.0 
time_units_string = "(-)"

# --- Output method for Paraview, HDF5 and  tracers 
# --- Every n steps or every time t ---
output_type = "steps" # "steps" or "time"
every_n 	= 1
every_t 	= 100*kyr

# --- Save tracers into files? ---
save_tracers = False

# --- What properties to save? Must be one of the following (order does not matter):
# rank              = to which process the tracer belongs
# dev_stress_xx     = xx component of the deviatoric stress
# dev_stress_xz     = xz component of the deviatoric stress
# plastic_strain    = amount of plastic strain
# ocean_material    = 1.0 if the tracer was marked as "ocean material"
# surface_material  = 1.0 if the tracer was marked as "surface material"
# original_depth    = original y-coordinate of the tracer
# composition_0     \
# composition_1       = material composition of the tracer
# composition_2     /
# melt_fraction     = amount of partial melt on the tracer
# origin            = 0 if the tracer is original, 1 if added later
# id                = unique ID of the tracer

Tracers_Output = ["composition_1", "rank", "composition_0"]

# --- Headers for the columns in the text file for tracers
# --- Up to the user (order corredponding to "stat_output").

Tracers_header = ["rank", "composition"]

# --- What functions to write into Paraview and HDF5 file ---
Paraview_Output = ["velocity", "composition_0", "tracers", "ranks"]

# --- What values to print in a text file every time step? Must be one of the following (order does not matter):
# nusselt 	= Nusselt number
# vrms 		= Root mean square velocity
# tracers 	= Number of tracers
# q_top 	= Heat flux over the top boundary
# q_bot 	= Heat flux over the bottom boundary
# time		= Duration of the simulation (hours)
# timestep	= Duration of the time step (seconds)

stat_output = ["vrms", "time", "timestep"]

# --- Headers for the columns in the text file
# --- Up to the user (order corredponding to "stat_output").

stat_header = ["v_rms (-)", "Time (h)", "dt (s)"]

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 2/ RELOADING OPTIONS -----------------------
#----------------------------------------------------------------------
reload_name 		= "van_Keken"

reloading_HDF5 		= False
reload_HDF5			= 12 # First column in data_timestamp.dat

reloading_tracers 	= False
reload_tracers		= 120 # Second column in data_timestamp.dat

# --- What functions to read from HDF5 file when reloading ---
Functions_Input = ["temperature"]

# --- Time step to reload (first column in data_timestamp.dat ---
restart_time 	= False

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------ 3/ MESH GEOMETRY ----------------------------
#----------------------------------------------------------------------
height = 1.0 # m
length = 0.9142 # m

loading_mesh = False
mesh_name 	= "meshes/mesh_25x100km.xml"

# --- Basic mesh resolution if not loading mesh ---
z_div = 50
x_div = int(z_div*(length/height)) # (keeps aspect ratio 1)
triangle_types = "crossed" # crossed, left, right, left/right, right/left

# --- Rectangular mesh refinement ---
# --- From x_left to x_right and y_bottom to y_top
# e.g. refinement = [x_left, x_right, y_bottom, y_top] ---

# --- Repeat within the [...] for multiple levels of refinement,
# leave empty for no refinement ---
refinement = []

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 4/ MATERIAL COMPOSITIOIN -------------------
#----------------------------------------------------------------------

# --- Simple material assignment using "interface", "circle" or "rectangle" geometries ---
# --- The n-th material overrides the (n-1)th.
# --- For more complicated geometries, modify "introduce_tracers" function in the m_tracers.py 

# --- Syntax: [[material_0], [material_1], ...] 
# Circle:    ["circle",     center_x,   center_y,   radius]
# Rectangle: ["rectangle",  left,       right,      bottom,     top]
# Interface: ["interface",  "below/above"]

def interface(x):
    return 0.02*cos(np.pi*x/length) + 0.2

# --- Leave empty for a single material ---
materials = [["interface", "below"], ["interface", "above"]]

# --- If True, the cells without tracers will be assigned material "default_composition" ---
# --- Applicable only if the material composition is the only tracer-requiring feature ---
allow_empty_cells = False
default_composition = 1
empty_region = [["rectangle",  0, 1, 0.3, 1]] #  So far rectangle implemented

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------ 4/ STOKES PROBLEM BOUNDARY CONDITIONS -------------------
#----------------------------------------------------------------------

# Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
BC_vel_top 		= "no_slip" 
BC_vel_bot 		= "no_slip"
BC_vel_left 	= "free_slip"
BC_vel_right 	= "free_slip"

# Influx boundary conditions - active only when "velocity" /_x /_y is selected above
# Top boundary
velocity_top_x = Constant(0.0)
velocity_top_y = Constant(0.0)

# Bottom boundary
velocity_bot_x = Constant(0.0)
velocity_bot_y = Constant(0.0)

# Left boudnary
velocity_left_x = Constant(-2.0e-11)
velocity_left_y = Constant(0.0)

# Right boundary
velocity_right_x = Constant(+2.0e-11)
velocity_right_y = Constant(0.0)

# Method for correcting velocity field in case of multiple free surface conditions
stokes_null = "volume"

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#------------ 5/ HEAT TRANSFER BOUNDARY CONDITIONS -------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# --- Solve heat transfer? ---
solve_energy_problem = True

BC_T_top 	= "temperature" 	# temperature, heat_flux or radiation
BC_T_bot 	= "temperature"		# temperature or heat_flux
BC_T_left 	= "heat_flux"		# temperature or heat_flux
BC_T_right 	= "heat_flux"		# temperature or heat_flux

# --- Boundary temperatures ---
T_top 	= 0.0	# K
T_bot 	= 1.0	# K
T_left 	= 0.0	# K
T_right = 0.0	# K

# --- Reference temperature ---
T_ref 	= 270.0 # K

# --- Heat fluxes in the direction of the boundary normal vector ---
bc_q_top 	= 0.0e-3	# W m^-2
bc_q_bot 	= 0.0e-3	# W m^-2
bc_q_left 	= 0.0e-3	# W m^-2
bc_q_right 	= 0.0e-3	# W m^-2

# --- Insolation parametes ---
emis 	= 0.97         	# Ice emissivity
albedo 	= 0.67
dist_AU = 5.2			# Object's distance from Sun in AU
insolation = insol_1AU/dist_AU**2

# --- Melting inside the shell ---
internal_melting = False
T_melt = 270.0	# K

# --- Tidal heating ---
tidal_dissipation = False
heating_model = "Maxwell" 	# Maxwell, Andrade or none
H_max = 2.717e-6			# W m^{-3}

# Andrade parameters
alpha_and = 0.2

# --- Find conductive initial condition ---
init_cond_profile = True

# Include tidal heating in the initial temperature?
initial_tidal_dissipation = False

# --- Cosine perturbation of initial temperature ---
cos_perturbation = True
perturb_ampl 	= 0.01 # K
perturb_freq 	= 0.5 # Number of cos waves in lateral direction
#----------------------------------------------------------------------

# === OUTPUT OPTIONS ===
t_end = 0.3 #20*Myr # s 



# --- Time scaling ---
# / deformation rate
# ----->
# x end time
# x constant timestep
# x topography diff. coef.
# x every_t

time_scaling = 1

# --- For surface and ocean material distribution ---
sm_thickness = 500.0	# m
om_thickness = 500.0	# m

tracers_per_cell = 25 # per reference area (value 25 is optimal)

weight_tracers_by_ratio = False

# --- Timestep strategy (domain, cell or  constant) ---
timestep_strategy = "domain"

# "domain" computes dt_cond based on (T_top + T_bot/2) 
# and dt_conv based on dx_min and v_max in the domain
dT_max = 4.0

# "cell" computes timestep in each cell and chooses the minimal

# "constant" prescribes a constant timestep
dt_const = 50*kyr #0.5*kyr*time_scaling


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#---------------------------- 6/ RHEOLOGY -----------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# === PHYSICAL INGREDIENTS ===
# --- Rheology ---
elasticity = False
plasticity = False

# --- Multiple composition ---
multiple_composition = False

# --- Phase transition at the bottom boundary ----
phase_transition = False
DAL_factor 	= 0.0 		# W/m3, strength of the phase transition
Lt 			= 334.0e3	# J/kg, latent heat

# Adjust oceanic heat flux in every time step?
keep_init_qw = False

# --- Adaptive topography diffusion and topography diffusion factor ---
adaptive_smoothing = False
topo_diff = 1e-8*time_scaling

# === PHYSICAL PARAMETERS ===
g 		= 1e4				# m s^-2
rho_s 	= 1.0 #920.0 #1.0 	# kg m^-3
rho_l 	= 1000.0 #2*Ra/(2.5e-5*1e3) - Ra*1.0		# kg m^-3
rho_m 	= rho_l		# kg m^-3
alpha_exp = 1.6e-4	# K^-1



# === STOKES PROBLEM ===
stokes_elements 	= "Mini" # Mini or TH (Taylor-Hood)
Picard_iter_max 	= 25
Picard_iter_error 	= 1e-3

# --- Advection parameters ---
cfl = 0.5
integration_method = "RK4" # Euler, RK2, RK4

# --- Rheological parameters ---
# --- VISCOSITY ---
viscosity_type = "GK_2001" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition

# --- Parameters for constant / temperature-dependent viscosity ---
eta_0 = 1.5e14	# Pa.s
Q_activ = 50e3 	# J/mol

# --- Parameters for nonlinear viscosity ---
d_grain = 1.0e-3

eta_max = 1e23
eta_min_plast = 1e13 #eta_max/1e6

# --- Elasticity ---
# shear_modulus = 7.0e10

# If elasticity is off, we can compute the new viscosity directly from the strain rate (it will be faster)
# However, the stress formula + VEP iterations would work too
viscosity_from = "strain_rate" # stress or strain_rate
stress_iter_error = 1e-4

# --- Plasticity ---
angle_matrix = 16.0*(3.1415926/180.0)
angle_inclusion = 0.0*(3.1415926/180.0)

C0 = 1e6
C_inf = 0.0

a0 = 16.0*(3.1415926/180.0)
a_inf = 4.57*(3.1415926/180.0)

eps_0 = 0
eps_inf= 0.2

healing = False
recovery_time = 300*kyr