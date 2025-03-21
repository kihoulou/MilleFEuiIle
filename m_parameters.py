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
# 3/ Simulation duration

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
name = "shear_bands"
"""
Name of the directory with the results. The directory with the results will be named ``data_name``.

:meta hide-value:
"""
override_directory = True # Protection from overwriting the directory above
"""
* **True**: Repeared run with the same ``name`` overwrites the data
* **False**: Original ``name`` directory will be protected, attempted run creates a directory ``data_name_new``

:meta hide-value:
"""
# --- In what units the time will be? ---
# 1.0 - seconds or nondimensional, or yr, kyr or Myr
time_units = 1.0
"""Time units at which the output will be written.

:var: Use ``1.0`` for seconds or in case of non-dimensional problems, or ``yr``, ``kyr`` or ``Myr``

:vartype: float

:meta hide-value:
"""
time_units_string = " "
"""
:var: String representation of ``time_units`` up to the user. 

:vartype: string

:meta hide-value:
"""

# --- Output method for Paraview, HDF5 and  tracers 
output_frequency = ["steps", 1] # ["steps", 10] or ["time", 100*kyr]
""" Frequency for the output to ParaView, HDF5 data and tracer files. The first time step is always saved, the
last one if the simulation arrives at :math:`t_{end}` succesfully.

:var:  
   * ``["steps", 10]`` outputs the data every 10 steps
   * ``["time", 100*kyr]`` outputs the data every 100 kyr

:vartype: list consisting of a string ``"steps"`` or ``"time"``, and an integer/float

:meta hide-value:
"""
# --- Save tracers into files? ---
save_tracers = False
"""
Whether to save tracers into files located at ``data_name/tracers``. Timing of the output is given by...

:meta hide-value:
"""
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
"""
List of strings indicating which properties carried by the tracers will be saved, e.g.::

   Tracers_Output = ["rank", "dev_stress_xx", "melt_fraction"]

.. table:: Properties to save on tracers
   :widths: auto

   =================  ==============
     String           Description
   =================  ==============
   rank                the process the tracer belongs to
   dev_stress_xx       *xx* component of the deviatoric stress
   dev_stress_xz       *xz* component of the deviatoric stress
   plastic_strain      amount of plastic strain
   ocean_material      1 if marked as "ocean material", 0 otherwise
   surface_material    1 if marked as "surface material", 0 otherwise
   original_depth      original *z*-coordinate of the tracer
   composition_n       concentration of material *n* (*n* being a number)
   melt_fraction       amount of partial melt on the tracer
   origin              0 if the tracer is original, 1 if added later
   id                  unique ID of the tracer
   =================  ==============

:meta hide-value:
"""
# --- Headers for the columns in the text file for tracers
# --- Up to the user (order corredponding to "stat_output").

Tracers_header = ["rank", "composition"]
"""
List of headers for the properties specified in ``tracers_output``, e.g.::

   Tracers_Output = ["rank (-)", "sigma_xx (Pa)", "xm (-)"]

:meta hide-value:
"""

# --- What functions to write into Paraview and HDF5 file ---
Paraview_Output = ["velocity",\
                "pressure",\
                "composition_0",\
                "composition_1",\
                "topography_top",\
                "strain_rate_inv",\
                "iteration_error",\
                "plastic_strain",\
                "stress_dev_inv",\
                "yield_stress",\
                "yield_function",\
                "eta_v",\
                "viscosity"]

# --- What values to print in a text file every time step? Must be one of the following (order does not matter):
# nusselt 	= Nusselt number
# vrms 		= Root mean square velocity
# tracers 	= Number of tracers
# q_top 	= Heat flux over the top boundary
# q_bot 	= Heat flux over the bottom boundary
# time		= Duration of the simulation (hours)
# timestep	= Duration of the time step (seconds)

stat_output = ["vrms", "time", "timestep", "avg_h_bot"]

# --- Headers for the columns in the text file
# --- Up to the user (order corredponding to "stat_output").

stat_header = ["v_rms (-)", "Time (h)", "dt (s)"]

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 2/ RELOADING OPTIONS -----------------------
#----------------------------------------------------------------------
reload_name 		= "van_Keken"
"""
:var: The directory from which the data will be reloaded when ``reloading_HDF5 = True``.

:vartype: string

:meta hide-value:
"""

reload_HDF5 		= False
"""
:var: Whether or not to load HDF5 data from ``data_reload_name/HDF5/data.h5``. 

:vartype: boolean

:meta hide-value:
"""

reload_HDF5_step			= 12 # First column in data_timestamp.dat
"""
:var: Time stamp of the HDF5 file which will be loaded when ``reloading_HDF5 = True``.

:vartype: integer

:meta hide-value:
"""

# --- What functions to read from HDF5 file when reloading ---
reload_HDF5_functions = ["temperature"]
"""
:var: Fields to be loaded from the HDF5 file ``data_reload_name/HDF5/data.h5``.

:vartype: list of strings

:meta hide-value:
"""

reload_tracers 	= False
"""
:var:  Whether to load tracers from ``data_reload_name/tracers/``.

:vartype: boolean

:meta hide-value:
"""

reload_tracers_step		= 120 # Second column in data_timestamp.dat
"""
:var:  Specifies the file from which to read the tracers ``data_reload_name/tracers/step_``.

:vartype: integer

:meta hide-value:
"""

# --- Time step to reload (first column in data_timestamp.dat ---
restart_time 	= False
"""
Whether to reset time when ``reloading_HDF5 = True``.

:var: 
   * **True** sets time to 0.
   * **False** continues from the given checkpoint.

:vartype: boolean

:meta hide-value:
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 3/ SIMULATION DURATION ---------------------
#----------------------------------------------------------------------
t_end = 0.02 #200*Myr # s 

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 3/ TRACERS OPTIONS -------------------------
#----------------------------------------------------------------------
sm_thickness = 500.0	# m
om_thickness = 500.0	# m

integration_method = "RK4" # Euler, RK2, RK4

tracers_per_cell = 20 

weight_tracers_by_ratio = False

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------ 3/ MESH GEOMETRY ----------------------------
#----------------------------------------------------------------------
height = 1 # m
length = 4 # m

loading_mesh = False
"""
:var: 
   * **True** reads a mesh from file ``mesh_name``.
   * **False** creates a new mesh based on ``x_div``, ``z_div`` and ``triangle_types`` 

:vartype: boolean

:meta hide-value:
"""

mesh_name 	= "meshes/mesh_25x100km.xml"

# --- Basic mesh resolution if not loading mesh ---
z_div = 100
x_div = int(z_div*(length/height)) # (keeps aspect ratio 1)
triangle_types = "crossed" # crossed, left, right, left/right, right/left
"""
Method of dividing basic squares into mesh triangle elements. 

:var: ``"crossed"``, ``"left"``, ``"right"``, ``"right/left"`` or ``"left/right"``
:vartype: string

:meta hide-value:
"""

# --- Rectangular mesh refinement ---
# --- From x_left to x_right and y_bottom to y_top
# e.g. refinement = [x_left, x_right, y_bottom, y_top] ---

# --- Repeat within the [...] for multiple levels of refinement,
# leave empty for no refinement ---
refinement = [] #[0.0, length, height/2.0, height]
""" 
:var: Minimum and maximum *x*- and *y*-coordinates of a rectangle area of the mesh to be refined::

      refinement = [0.0, length, height/2.0, height]

   Repeating the sequence superposes the refinements. For example::

      refinement = [0.0, length, 0, height/2.0, 0.0, length, 0, height/4.0]

   gives first-level refinement in the bottom half domain and second-level one in the bottom quarter of the domain.

:vartype: list of 4x *n* floats, where *n* is the number of refined areas
   
:meta hide-value:
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 4/ MATERIAL COMPOSITIOIN -------------------
#----------------------------------------------------------------------

# --- Simple material assignment using "interface", "circle" or "rectangle" geometries ---
# --- The n-th material overrides the (n-1)th.
# --- For more complicated geometries, modify "introduce_tracers" function in the m_tracers.py 

# --- Syntax: [[material_0], [material_1], ...] 
# Circle:               ["circle",     center_x,   center_y,   radius]
# Rectangle:            ["rectangle",  left,       right,      bottom,     top]
# Interface:            ["interface",  "below/above"]

def interface(x):
    """A function just for me.

    .. note:: This is a note admonition.
    .. math:: f(x) = Acos(x)


    :param my_arg: The first of my arguments.
    :param my_other_arg: The second of my arguments.

    :returns: A message (just for me, of course).
    """
    # return 0.02*cos(np.pi*x/length) + 0.2
    return 5e2*cos(2.0*np.pi*x/length) + 3e3

# --- Leave empty for a single material ---
dd = 0.08
ll = length/2.0 - dd/2.0
rr = length/2.0 + dd/2.0
bb = 0.0
tt = dd/2.0
materials = [["rectangle",  0.0, length, 0.0, height], ["rectangle",  ll, rr, bb, tt]] # [["interface", "above"], ["interface", "below"]]

# --- If True, the cells without tracers will be assigned material "default_composition" ---
# --- Applicable only if the material composition is the only tracer-requiring feature ---
allow_empty_cells = False
"""
:var: Whether or not allow to have elements without tracers
:vartype: boolean

:meta hide-value:
"""

default_composition = 0
"""
:var: If ``allow_empty_cells = True`` this determines the composition of empty cells
:vartype: integer

:meta hide-value:
"""

empty_region = [["rectangle",  0, length, 20e3, height]] #  So far rectangle implemented
"""
:var: Specifies geometry of the region where the empty cells will be be given  ``default composition`` even without tracers inside.
:vartype: integer

:meta hide-value:
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------ 4/ STOKES PROBLEM BOUNDARY CONDITIONS -------------------
#----------------------------------------------------------------------
stokes_elements = "TH"

time_step_position = "stokes" #stokes (right after Stokes problem) or "end" (at the end of the time loop)
"""Determines where to compute new time step.

:var: ``"stokes"`` after the Stokes solver or  ``"end"`` at the end of the time loop
:vartype: string

:meta hide-value:
"""

time_step_strategy = "constant"
""" Method for computing a new time step

:var:
   
   * ``"domain"`` computes a new time step based on specific criteria in every step of the time loop (see ``m_timestep`` module)
   * ``"constant"`` prescribes a constant time step given by xxx

:vartype: string

:meta hide-value:
"""

cfl = 0.5
Picard_iter_max 	= 25
""" 
:var: Maximum number of Stokes solver Picard iterations in case of nonlinear rheology (stress-dependent viscosity, plasticity).
:vartype: integer

:meta hide-value:
"""

Picard_iter_error 	= 1e-3
""" 
:var: Minimum relative error in velocity field sufficient to exit the Stokes solver Picard iterations.
:vartype: float

:meta hide-value:
"""

error_type          = "maximum" # "maximum or integrated"

# Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
BC_vel_top 		= "free_surface" 
BC_vel_bot 		= "free_slip"
BC_vel_left 	= "velocity_x"
BC_vel_right 	= "velocity_x"

mesh_displacement_laplace = "full"  #"full_laplace" or "z_only"
""" Determines whether full Laplace equation or only its *z*-part will 
be solved to obtain mesh displacement on the mesh between the top and bottom boundary.

:var: 
   * ``"full"`` distributes the displacement by solving :math:`\\nabla^2 h_2 = 0`\ . Useful for simulation involving tectonics, mesh is more stable.
   * ``"z_only"`` distributes the displacement linearly by solving :math:`\\nabla^2_z h_2 = 0`

:vartype: string

:meta hide-value:
"""

# Influx boundary conditions - active only when "velocity" /_x /_y is selected above
# Top boundary
velocity_top_x = Constant(0.0)
velocity_top_y = Constant(0.0)

# Bottom boundary
velocity_bot_x = Constant(0.0)
velocity_bot_y = Constant(0.0)

# Left boundary
velocity_left_x = Constant(+2.0)
velocity_left_y = Constant(0.0)

# Right boundary
velocity_right_x = Constant(-2.0)
velocity_right_y = Constant(0.0)

# Method for correcting velocity field in case of multiple free surface conditions
stokes_null = "boundary"

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#------------ 5/ HEAT TRANSFER BOUNDARY CONDITIONS -------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# --- Solve heat transfer? ---
solve_energy_problem = False
nonlinear_heat_equation = False

BC_T_top 	= "temperature" 	# temperature, heat_flux or radiation
BC_T_bot 	= "temperature"		# temperature or heat_flux
BC_T_left 	= "heat_flux"		# temperature or heat_flux
BC_T_right 	= "heat_flux"		# temperature or heat_flux

# --- Boundary temperatures ---
T_top 	= 90.0	# K
T_bot 	= 270.0	# K
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
perturb_ampl 	= 1 # K
perturb_freq 	= 1 # Number of cos waves in lateral direction
#----------------------------------------------------------------------

# --- Time scaling ---
# / deformation rate
# ----->
# x end time
# x constant timestep
# x topography diff. coef.
# x every_t

time_scaling = 1

# "domain" computes dt_cond based on (T_top + T_bot/2) 
# and dt_conv based on dx_min and v_max in the domain
dT_max = 4.0

# "cell" computes timestep in each cell and chooses the minimal

# "constant" prescribes a constant timestep
dt_const = 0.0005 #0.5*kyr #0.5*kyr*time_scaling

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#---------------------------- 6/ RHEOLOGY -----------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# === PHYSICAL PARAMETERS ===
g 		= 1.0				# m s^-2
rho_s 	= 2700.0 #1.0 	# kg m^-3
rho_l 	= 1000.0 #2*Ra/(2.5e-5*1e3) - Ra*1.0		# kg m^-3
rho_m 	= rho_l		# kg m^-3
alpha_exp = 1.6e-4	# K^-1

# === PHYSICAL INGREDIENTS ===
# --- Rheology ---
elasticity = False
plasticity = True

# --- Phase transition at the bottom boundary ----
phase_transition = False
DAL_factor 	= 1e-2		# W/m3, strength of the phase transition
Lt 			= 334.0e3	# J/kg, latent heat

# --- Adaptive topography diffusion and topography diffusion factor ---
adaptive_smoothing = False
topo_diff = 1e-8*time_scaling

# --- Rheological parameters ---
# --- VISCOSITY ---
viscosity_type = "composition" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition

# --- Parameters for constant / temperature-dependent viscosity ---
eta_0 = 1e15	# Pa.s
Q_activ = 60e3 	# J/mol

# --- Parameters for nonlinear viscosity ---
d_grain = 1.0e-3

eta_max = 1e23
eta_min_plast = 1e-13 #eta_max/1e6

# --- Elasticity ---
shear_modulus = 3.52e9

# If elasticity is off, we can compute the new viscosity directly from the strain rate (it will be faster)
# However, the stress formula + VEP iterations would work too
stress_iter_error = 1e-4

# --- Plasticity ---
int_friction_angle = 40.0*(np.pi/180.0)
int_friction_angle2 = 0.0*(np.pi/180.0)
"""
Angle of internal friction in radians.

:meta hide-value:
"""

cohesion_strong, cohesion_weak = 400.0, 20.0
"""
Cohesion of an undamaged and a fully damaged material, respectively.

:meta hide-value:
"""

yield_stress_max, yield_stress_min = 1e4, 200
"""
Upper and lower cut-off values for the yield stress, respectively.

:meta hide-value:
"""

eps_strong = 0
"""
Value of the plastic strain at which the material starts to accumulate damage.

:meta hide-value:
"""

eps_weak = 0.1
"""
Critical value of the plastic strain beyond which the material is considered as fully damaged.

:meta hide-value:
"""

healing = False
"""
Whether the accumulated plastic strain decreases in time due to microscopic healing processes.
If ``True`` the the plastic strain is decreased by a following formula

.. math::

:meta hide-value:
"""

recovery_time = 300*kyr
"""
Characteristic time scale for the microscopic healing processes.

:meta hide-value:
"""