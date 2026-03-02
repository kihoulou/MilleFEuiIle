from dolfin import *
from m_constants import *
import numpy as np
import sys

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------- 1/ OUTPUT FILES SETTINGS -------------------
#----------------------------------------------------------------------
# --- Name of the directory with results ---
name = "demo_convection"

# Protection from overwriting the directory above
protect_directory = False

# --- Name of the directory from which the data will be reloaded ---
reload_name 		= ""

# --- Whether or not to load HDF5 data from reload_name/HDF5/data.h5 --- 
reload 		= False

# ---  Time stamp of the HDF5 file which will be loaded ---
reload_step			= 10

# --- Whether to reset time when reloading ---
restart_time 	= False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------------  TIME SETTINGS -------------------
#----------------------------------------------------------------------

# --- In what units the time will be? ---
# 1.0 - seconds or nondimensional, or yr, kyr or Myr
time_units = kyr

# --- String representation of time_units ---
time_units_string = "kyr"

# --- Criterion for ending the simulation, e.g., ["time", 1*Myr], ["step", 1000] or ["initial_condition", *]---
termination_condition = ["time", 10*Myr]
# termination_condition = ["step", 30]

# --- Output method for Paraview, HDF5 and  tracers ---
output_frequency = ["steps", 10] #["time", 50*kyr] # e.g. ["steps", 10] or ["time", 100*kyr]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------- FILE SETTINGS -------------------
#----------------------------------------------------------------------

# --- What properties to save? Must be one of the following (order does not matter):
# rank              = to which process the tracer belongs
# dev_stress_xx     = xx component of the deviatoric stress
# dev_stress_xz     = xz component of the deviatoric stress
# plastic_strain    = amount of plastic strain
# original_depth    = original y-coordinate of the tracer
# composition_0     \
# composition_1       = material composition of the tracer
# composition_2     /
# melt_fraction     = amount of partial melt on the tracer
# origin            = 0 if the tracer is original, 1 if added later
# id                = unique ID of the tracer
Tracers_Output = []

# --- Headers for the columns in the text file for tracers ---
Tracers_header = [] # KEEP EMPTY, will be filled automatically

# --- What functions to write into Paraview and HDF5 file ---
Paraview_Output = ["temperature", "viscosity", "velocity", "iteration_error", "pressure"]

Paraview_Output_Ini = ["temperature"]

# --- What values to print in a text file every time step? Must be one of the following (order does not matter):
# nusselt 	= Nusselt number
# vrms 		= Root mean square velocity
# tracers 	= Number of tracers
# q_top 	= Heat flux over the top boundary
# q_bot 	= Heat flux over the bottom boundary
# time		= Duration of the simulation (hours)
# timestep	= Duration of the time step (seconds)
stat_output = ["q_top", "time", "timestep"]

# --- Headers for the columns in the text file
# --- Up to the user (order corredponding to "stat_output").
stat_header = ["q_top (mW/m2)", "Time (h)", "dt (s)"]

monitor_cache = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------ 3/ MESH SETTINGS ----------------------------
#----------------------------------------------------------------------

# --- Mesh height ---
height = 167e3 # m

# --- Mesh length ---
length = 167e3*2 # m

# --- Whether to read an external mesh ---
loading_mesh = False

# --- Mesh to be loaded ---
mesh_name 	= ""

# --- Basic mesh resolution if not loading mesh ---
z_div = 50

# --- Number of nodes in horizontal direction ---
x_div = int(z_div*(length/height)) # (keeps aspect ratio 1)

# --- Method of dividing basic squares into mesh triangle elements. ---
triangle_types = "crossed" # crossed, left, right, left/right, right/left

# --- Rectangular mesh refinement ---
# --- From x_left to x_right and y_bottom to y_top
# e.g. refinement = [x_left, x_right, y_bottom, y_top] ---

# --- Repeat within the [...] for multiple levels of refinement,
# leave empty for no refinement ---
refinement = [] #[length/4.0, 3*length/4.0, 0, height]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------- 4/ MATERIAL COMPOSITION -------------------
#----------------------------------------------------------------------

# --- Simple material assignment using "interface", "circle" or "rectangle" geometries ---
# --- The n-th material overrides the (n-1)th.
# --- For more complicated geometries, modify "introduce_tracers" function in the m_tracers.py 

# --- Syntax: [[material_0], [material_1], ...] 
# Circle:               ["circle",     center_x,   center_y,   radius]
# Rectangle:            ["rectangle",  left,       right,      bottom,     top]
# Interface:            ["interface",  "below/above"]

def interface(x):
    return 75e3*exp(-(x - 100e3)**2/1e8)

# --- Leave empty for a single material ---
materials = []

# --- If True, the cells without tracers will be assigned material "default_composition" ---
# --- Applicable only if the material composition is the only tracer-requiring feature ---
empty_cells_allowed = False

empty_cells_composition = 0

empty_cells_region = [] #  So far rectangle implemented

initial_topography = False

h_top_ini = Expression("0.0*sin(120*pi*x[0]/l)", l = length, pi = np.pi, degree = 1)

h_bot_ini = Expression("8.5e3*(tanh((x[0]-l/4)/3000) - tanh((x[0] - 3*l/4)/3000))/2.0", l = length, degree = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------ 4/ STOKES PROBLEM SETTINGS -------------------
#----------------------------------------------------------------------
stokes_elements = "Mini"

time_step_position = "stokes"# (right after Stokes problem) or "end" (at the end of the time loop)

time_step_strategy = "domain"

cfl = 0.5

# ---  Maximum number of Stokes solver Picard iterations (used for heat equation, advection of tracers) ---
Picard_iter_max 	= 10

# --- Minimum relative error in velocity field ---
Picard_iter_error 	= 1e-3

error_type          = "integrated" # "maximum or integrated"

# Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
BC_Stokes_problem = [["free_slip"], # top boundary (1)
                     ["free_slip"], # bottom boundary (2)
                     ["free_slip"], # left boundary (3)
                     ["free_slip"]] # right boundary (3)

mesh_displacement_laplace = "full"  #"full" or "z_only"

# Method for correcting velocity field in case of multiple free surface conditions
stokes_null = "boundary"

integration_method = "RK4" # Euler, RK2, RK4

tracers_per_cell = 15

weight_tracers_by_ratio = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#:::::::::::::::::::: HEAT TRANSFER SETTINGS ::::::::::::::::::::::::::
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Solve heat transfer? ---
solve_energy_problem = True

# --- Whether the heat transfer should be solved with nonlinear solver ---
nonlinear_heat_equation = False

# --- Boundary condition for heat transfer equation ---
BC_heat_transfer   = [["temp", 100.0],     # top boundary    (1)
                     ["temp", 270.0],     # bottom boundary (2)
                     ["heat_flux", 0.0],  # left boundary   (3)
                     ["heat_flux", 0.0]]  # right boundary  (4)

# --- Reference temperature ---
T_ref 	= 270.0 # K

# --- Insolation parametes ---
emis 	= 0.97         	# Ice emissivity
albedo 	= 0.67
dist_AU = 5.2			# Object's distance from Sun in AU
insolation = insol_1AU/dist_AU**2

# --- Melting inside the shell ---
# --- Whether generate partial melt if the temperature of the solid reaches T_melt ---
internal_melting = False

# --- Melting temperature of the solid ---
T_melt = 270.0	# K

# --- Threshold for adding or removing tracers when melting ---
dT_melt_tresh = 4.0

# --- Tidal heating ---
tidal_dissipation           = False

initial_tidal_dissipation   = False
heating_model               = "Maxwell" # Maxwell, Andrade or none
H_max = 4e-6 # W m^{-3}

# Andrade parameters
alpha_and = 0.2

# --- Find conductive initial condition ---
init_cond_profile = True

# --- Cosine perturbation of initial temperature ---
cos_perturbation = True

# --- Thermal perturbation amplitude ---
perturb_ampl 	= 5 # K
# --Number of half-cosine waves in lateral direction ---
perturb_freq 	= 1 

# "domain" computes dt_cond based on (T_top + T_bot/2) 
# and dt_conv based on dx_min and v_max in the domain
dT_max = 4.0

# "cell" computes timestep in each cell and chooses the minimal

# "constant" prescribes a constant timestep
dt_const = kyr

maximum_time_step = False
time_step_scaling = 25 # dt = t / time_step_scaling

scaled_time_step = False
dt_max = 0.1*Myr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---------------------------- 6/ RHEOLOGY -----------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# === PHYSICAL PARAMETERS ===
# --- Density of the domain solid ---
rho_s 	= 920.0 # kg/m3

# --- Density of the liquid below the domain ---
rho_l 	= 1000.0 # kg/m3

# --- Density of the melt produced by the solid ---
rho_m 	= rho_l # kg/m3

alpha_exp = 1.6e-4	# K^-1

# === PHYSICAL INGREDIENTS ===
# --- Gravity acceleration ---
g 		= 1.3				# m s^-2

# --- Rheology ---
elasticity = False

G_ice = 3.52e9   # Pa

plasticity = False

# --- Phase transition at the bottom boundary ----
phase_transition = False

# --- Strength of the phase transition at the ice-water boundary ---
DAL_factor = 1e-3 # W/m3

# --- Latent heat of the material ---
Lt 			= 334.0e3 # J/kg

# --- Adaptive topography diffusion and topography diffusion factor ---
adaptive_smoothing = False
topo_diff = 1e-8

# --- Rheological parameters ---
# --- VISCOSITY ---
viscosity_type = "temp-dep" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition

# --- Parameters for constant / temperature-dependent viscosity ---
eta_0 = 3e15	# Pa.s

# --- Activation energy ---
Q_activ = 60e3 	# J/mol

# --- Grain size ---
d_grain = 1.0e-3

# --- Upper cut-off viscosity ---
eta_max = 1e24

# --- Lower cut-off value for plastic viscosity ---
eta_min_plast = 1e10 #eta_max/1e6
eta_min = 1e10 #eta_max/1e6

# # --- Elasticity ---
# If elasticity is off, we can compute the new viscosity directly from the strain rate (it will be faster)
# However, the stress formula + VEP iterations would work too
stress_iter_error = 1e-4

# --- Plasticity ---
# --- Angle of internal friction in degrees ---
int_friction_angle = 16.0
int_friction_angle2 = 0.0

# --- Cohesion of an undamaged material (Pa) ---
cohesion_strong = 1e6

# --- Cohesion of a fully damaged material ---
cohesion_weak = 0.0

yield_stress_max, yield_stress_min = 1e8, 0.1

# --- Value of the plastic strain at which the material starts to accumulate damage ---
eps_strong = 0

# --- Critical value of the plastic strain beyond which the material is considered as fully damaged ---
eps_weak = 0.2

# --- Whether the accumulated plastic strain decreases in time ---
healing = False

# --- Characteristic time scale for the microscopic healing processes ---
healing_timescale = 100*kyr
