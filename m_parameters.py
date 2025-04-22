from dolfin import *
from m_constants import *
import numpy as np
import sys

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
name = "rising_plume_" + str(sys.argv[1]) + "_eta" + str(sys.argv[2])
"""
:var: Name of the directory with the results. The directory with the results will be named ``data_name``.


:vartype: string

:meta hide-value:
"""

# Protection from overwriting the directory above
protect_directory = False
"""
* **True**: Repeated run protects the original data by creating saving into a directory ``data_name_new`` instead
* **False**: Repeated run with the same ``name`` overwrites the data

:meta hide-value:
"""
# --- In what units the time will be? ---
# 1.0 - seconds or nondimensional, or yr, kyr or Myr
time_units = Myr
"""The units of time at which the output will be written.

:var: Use ``1.0`` for seconds or in case of non-dimensional problems, or ``yr``, ``kyr`` or ``Myr``

:vartype: float

:meta hide-value:
"""

# --- String representation of time_units ---
time_units_string = "Myr"
"""
:var: String representation of ``time_units``, up to the user. 

:vartype: string, e.g.  ``"s"``, ``"yr"``, ``"kyr"`` or ``"Myr"``

:meta hide-value:
"""

# --- Output method for Paraview, HDF5 and  tracers ---
output_frequency = ["time", 1*Myr] # e.g. ["steps", 10] or ["time", 100*kyr]
""" 
:var: Frequency for the output to ParaView, HDF5 data and mesh, and tracer files. The first time step is always saved, the
      last one if the simulation succesfully ends through the termination criteria.

:vartype: list consisting of a string ``"steps"`` or ``"time"``, and an integer/float
   For example::

      output_frequency = ["steps", 10]

   outputs the data every 10 steps, while::

      output_frequency = ["time", 100*kyr]
   
   outputs the data every 100 kyr.

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
                "density",\
                "composition_0",\
                "composition_1",\
                "topography_top",\
                "strain_rate_inv",\
                "stress_dev_inv",\
                "tracers",\
                "z_function",\
                "shear_modulus",\
                "viscosity"]

# --- What values to print in a text file every time step? Must be one of the following (order does not matter):
# nusselt 	= Nusselt number
# vrms 		= Root mean square velocity
# tracers 	= Number of tracers
# q_top 	= Heat flux over the top boundary
# q_bot 	= Heat flux over the bottom boundary
# time		= Duration of the simulation (hours)
# timestep	= Duration of the time step (seconds)
stat_output = ["h_top_max", "time", "timestep"]

# --- Headers for the columns in the text file
# --- Up to the user (order corredponding to "stat_output").
stat_header = ["h_top (m)", "Time (h)", "dt (s)"]

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 2/ RELOADING OPTIONS -----------------------
#----------------------------------------------------------------------

# --- Name of the directory from which the data will be reloaded ---
reload_name 		= "van_Keken"
"""
:var: The directory from which the data will be reloaded when ``reloading_HDF5 = True``.

:vartype: string

:meta hide-value:
"""

# --- Whether or not to load HDF5 data from data_reload_name/HDF5/data.h5 --- 
reload_HDF5 		= False
"""
:var: Whether or not to load HDF5 data from ``data_reload_name/HDF5/data.h5``. 

:vartype: boolean

:meta hide-value:
"""

# ---  Time stamp of the HDF5 file which will be loaded ---
reload_HDF5_step			= 12
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

# --- Whether to reset time when reloading ---
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

# --- Criterion for ending the simulation, e.g. ["time", 1*Myr] or ["step", 1000] ---
termination_condition = ["time", 20*Myr]
""" Time or step criterion for ending the simulation naturally.

:var:  
   * ``["step", 1000]`` exits the time loop after the 100th step
   * ``["time", 100*kyr]`` exits the time loop as soon as 100 kyr is reached 


:meta hide-value:
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------- 3/ TRACERS OPTIONS -------------------------
#----------------------------------------------------------------------
sm_thickness = 500.0	# m
om_thickness = 500.0	# m

integration_method = "RK4" # Euler, RK2, RK4
""" Method for integration of the trajectories of Lagrangian tracers.

:var:  
   * ``"RK2"`` for the 2nd order Runge-Kutta scheme
   * ``"RK4"`` for the 4nd order Runge-Kutta scheme

:vartype: string

:meta hide-value:
"""

tracers_per_cell = 20 

weight_tracers_by_ratio = False

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------------------ 3/ MESH GEOMETRY ----------------------------
#----------------------------------------------------------------------

# --- Mesh height ---
height = 700e3 # m
""" Height of the rectangular mesh.

:vartype: float

:meta hide-value:
"""

# --- Mesh length ---
length = 2800e3 # m
""" Length of the rectangular mesh.

:vartype: float

:meta hide-value:
"""

# --- Whether to read an external mesh ---
loading_mesh = False
""" Whether to read an external mesh.
:var: 
   * **True** reads a mesh from file ``mesh_name``.
   * **False** creates a new mesh based on ``x_div``, ``z_div`` and ``triangle_types`` 

:vartype: boolean

:meta hide-value:
"""

# --- Mesh to be loaded ---
mesh_name 	= "meshes/mesh_25x100km.xml"
""" Mesh to be loaded.

:var: path to a mesh in **xml** format, e.g. ``mesh_25x100km.xml``

:vartype: string

:meta hide-value:
"""

# --- Basic mesh resolution if not loading mesh ---
z_div = 8
""" Number of nodes in vertical direction.

:vartype: integer

:meta hide-value:
"""

# --- Number of nodes in horizontal direction ---
x_div = int(z_div*(length/height)) # (keeps aspect ratio 1)
""" Number of nodes in horizontal direction.

:var: default ``int(z_div*(length/height))`` which keeps acpect ratio of the elements equal to 1

:vartype: integer

:meta hide-value:
"""

# --- Method of dividing basic squares into mesh triangle elements. ---
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
refinement = [0, length, height - 200e3, height,\
			  1200e3, 1600e3, 0, height-250e3,\
			  0, length, height - 200e3, height,\
			  1200e3, 1600e3, 0, height-220e3,\
			  0, length, height - 160e3, height,\
			  1250e3, 1550e3, 160e3, height-165e3]
# ,\
         #   1250e3, 1550e3, 200e3, height]
""" 
:var: Minimum and maximum *x*- and *y*-coordinates of a rectangle area of the mesh to be refined, see :func:`m_mesh.MeshModule.refine_mesh`\ .

:vartype: list of 4x *n* floats, where *n* is the number of refined areas

.. hint:: 
   In order to refine the upper half of the mesh: ::

      refinement = [0.0, length, height/2.0, height]
   
   By repeating the sequence, the areas of refinement can be superposed. For example ::

      refinement = [0.0, length, 0, height/2.0, 0.0, length, 0, height/4.0]

   gives first-level refinement in the bottom half domain and second-level one in the bottom quarter of the domain.

.. warning:: For tectonic simulations, the top boundary nodes need to be **equidistant** because of 
   detecting and smoothing node oscillations, see :func:`m_mesh.detect_oscillations`\ .

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
eta_mantle, eta_lid, eta_plume = 1e21, 10.0**int(sys.argv[2]), 1e20
G_mantle, G_lid, G_plume = 1e20, 7e10, 1e20
rho_mantle, rho_lid, rho_plume = 3300.0, 3300.0, 3200.0

materials = [["rectangle",  0.0, length, 0.0, height], ["rectangle",  0.0, length, 600e3, 700e3], ["circle", 1400e3, 300e3, 50e3]]

# --- If True, the cells without tracers will be assigned material "default_composition" ---
# --- Applicable only if the material composition is the only tracer-requiring feature ---
empty_cells_allowed = False
"""
:var: Whether to allow elements without tracers. If ``True``, tracers will not be added
      in :func:`m_tracers.Tracers.introduce_tracers` nor :func:`m_tracers.Tracers.add_tracers`\ .
:vartype: boolean

.. warning:: Use ``True`` only when **material composition** is the only tracer-requiring feature.

:meta hide-value:
"""

empty_cells_composition = 0
"""
:var: If ``empty_cells_allowed = True``, this determines the composition of empty cells,
      see :func:`m_interpolation.composition_interpolation`.
:vartype: integer

:meta hide-value:
"""

empty_cells_region = [["rectangle",  0, length, 20e3, height]] #  So far rectangle implemented
"""
:var: Specifies geometry of the region where the empty cells will be given  ``empty_cells_composition`` even without tracers inside,
      see :func:`m_tracers.Tracers.introduce_tracers`\ .
:vartype: list of lists consisting of the region shape and its dimensions.
   For example ::

      empty_cells_region = [["rectangle",  0, length, 20e3, height]]

.. note:: Only ``"rectangle"`` type is implemented at this point.

:meta hide-value:
"""

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#------------ 4/ STOKES PROBLEM BOUNDARY CONDITIONS -------------------
#----------------------------------------------------------------------
stokes_elements = "Mini"

time_step_position = "end" #stokes (right after Stokes problem) or "end" (at the end of the time loop)
"""Determines where to compute new time step.

:var: ``"stokes"`` after the Stokes solver or  ``"end"`` at the end of the time loop
:vartype: string

:meta hide-value:
"""

time_step_strategy = "domain"
""" Specifies the method for computing a new time step, see :func:`m_timestep.time_step`\ .

:var:
   
   * ``"domain"`` computes a new time step based on specific criteria in every step of the time loop (see ``m_timestep`` module)
   * ``"constant"`` prescribes a constant time step given by xxx

:vartype: string

.. tip:: For tectonic simulations, a time step smaller than given by the CFL criterion needs to be used to ensure sharp faults.

:meta hide-value:
"""

# --- The CFL parameter ---
cfl = 0.5
""" 
:var: The CFL parameter, determining the length of the adaptive time step, see :func:`m_timestep.time_step`\ .
:vartype: float

:meta hide-value:
"""

# ---  Maximum number of Stokes solver Picard iterations ---
Picard_iter_max 	= 25
""" 
:var: Maximum number of Stokes solver Picard iterations in case of nonlinear rheology (stress-dependent viscosity, plasticity),
      see :func:`m_equations.Equations.solve_Stokes_problem`\ .
:vartype: integer

:meta hide-value:
"""

# --- Minimum relative error in velocity field ---
Picard_iter_error 	= 1e-3
""" 
:var: Minimum relative error in velocity field sufficient to exit the Stokes solver Picard iterations,
      see :func:`m_equations.Equations.solve_Stokes_problem`\ .
:vartype: float

:meta hide-value:
"""

error_type          = "maximum" # "maximum or integrated"

# Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
BC_Stokes_problem = [["free_surface"],                # top boundary       (1)
                     ["no_slip"],    # bottom boundary    (2)
                     ["free_slip"],       # left boundary      (3)
                     ["free_slip"]]       # right boundary     (4)
"""
:var: Boundary conditions for the Stokes problem :func:`m_equations.Equations.equation_Stokes` in the following order

   * top boundary 
      * ``["free_slip/no_slip/free surface"]``,
      * ``["velocity_x/velocity_y", value]`` for specifying one component of the velocity,
      * or ``["velocity", value_x, value_y]`` for specifying both components of the velocity
   * bottom boundary 
      * ``["free_slip/no_slip/free surface"]``,
      * ``["velocity_x/velocity_y", value]`` for specifying one component of the velocity,
      * or ``["velocity", value_x, value_y]`` for specifying both components of the velocity
   * left boundary
      * ``["free_slip/no_slip"]``,
      * ``["velocity_x/velocity_y", value]`` for specifying one component of the velocity,
      *  or ``["velocity", value_x, value_y]`` for specifying both components of the velocity
   * right boundary
      * ``["free_slip/no_slip"]``,
      * ``["velocity_x/velocity_y", value]`` for specifying one component of the velocity,
      * or ``["velocity", value_x, value_y]`` for specifying both components of the velocity

:vartype: 
   * list of four lists (one for each boundary)
   * unit m s\ :math:`^{-1}` 

.. hint:: 
   For example, for a tectonic simulation with free surface at the top boundary, ouflow through the side boundaries by 10 km/Myr
   and strictly vertical inflow through the bottom boundary at the same rate, the boundary condition will look as follows: ::

      BC_Stokes_problem = [["free_surface"],          # top boundary       (1)
                     ["velocity", 0.0, +10e3/Myr],    # bottom boundary    (2)
                     ["velocity_x", -10e3/Myr],       # left boundary      (3)
                     ["velocity_x", +10e3/Myr]]       # right boundary     (4)

:meta hide-value:
"""

mesh_displacement_laplace = "full"  #"full_laplace" or "z_only"
""" Determines whether full Laplace equation or only its *z*-part will 
be solved to obtain mesh displacement on the mesh between the top and bottom boundary.

:var: 
   * ``"full"`` distributes the displacement by solving :math:`\\nabla^2 h_2 = 0`\ . Useful for simulation involving tectonics, mesh is more stable.
   * ``"z_only"`` distributes the displacement linearly by solving :math:`\\nabla^2_z h_2 = 0`

:vartype: string

:meta hide-value:
"""

# Method for correcting velocity field in case of multiple free surface conditions
stokes_null = "boundary"
""" 
:var: Method for correcting velocity field in case of free surface conditions on both top and bottom boundary.
:vartype:
   * ``"boundary"`` for subtracting the average :math:`v_z`\ integrated along the top boundary 
   * ``"volume"`` for subtracting the average :math:`v_z`\ integrated over the full domain

:meta hide-value:
"""

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#------------ 5/ HEAT TRANSFER BOUNDARY CONDITIONS -------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# --- Solve heat transfer? ---
solve_energy_problem = False
"""
:var: Whether to solve heat transfer :func:`m_equations.Equations.equation_heat`\ .

:vartype: boolean

:meta hide-value:
"""

# --- Whether the heat transfer should be solved with nonlinear solver ---
nonlinear_heat_equation = False
"""
:var: Whether the heat transfer :func:`m_equations.Equations.equation_heat` should be solved with nonlinear solver.

:vartype: boolean

:meta hide-value:
"""

# --- Boundary condition for heat transfer equation ---
BC_heat_transfer   = [["temp", 90.0],     # top boundary    (1)
                     ["temp", 270.0],     # bottom boundary (2)
                     ["heat_flux", 0.0],  # left boundary   (3)
                     ["heat_flux", 0.0]]  # right boundary  (4)
"""
:var: Boundary conditions for heat transfer equation in the following order

   * top boundary ``["temp/heat_flux/radiation", value]`` 
   .. note:: If ``"radiation"`` boundary condition is selected, the ``value`` serves as an initial estimate of the top boudnary temperature for the nonlinear solver.
   * bottom boundary ``["temp/heat_flux", value]``
   * left boundary ``["temp/heat_flux", value]``
   * right boundary ``["temp/heat_flux", value]`` 

:vartype: list of four lists (one for each boundary)

   * unit K for temperature
   * unit W m\ :math:`^{-2}` for the heat flux

.. hint:: 
   For example, in order to prescribe temperature 90 K at the top boundary, 265 K at the bottom boundary and insulating side boundaries,
   the boundary condition will look as follows: ::

      BC_heat_transfer = [["temp", 90.0],    # top boundary    (1)
                        ["temp", 265.0],     # bottom boundary (2)
                        ["heat_flux", 0.0],  # left boundary   (3)
                        ["heat_flux", 0.0]]  # right boundary  (4)

:meta hide-value:
"""

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
"""
:var: Whether generate partial melt if the temperature of the solid reaches :math:`T_{melt}`\ .

:vartype: boolean

:meta hide-value:
"""

# --- Melting temperature of the solid ---
T_melt = 270.0	# K
"""
:var: Melting temperature of the solid.

:vartype: float, unit K

:meta hide-value:
"""

# --- Tidal heating ---
tidal_dissipation           = False
"""
:var: Whether there is a volumetric source term in the heat transfer equation.

:vartype: boolean
:meta hide-value:
"""

initial_tidal_dissipation   = False
heating_model               = "Maxwell" # Maxwell, Andrade or none
H_max = 4e-6 # W m^{-3}

# Andrade parameters
alpha_and = 0.2

# --- Find conductive initial condition ---
init_cond_profile = False
"""
:var: Whether to find an initial conductive profile as an initial condition for the heat transfer problem.

:vartype: boolean

.. note:: If ``False``, then a linear temperature profile between the top and bottom temperature will be prescribed.

:meta hide-value:
"""

# --- Cosine perturbation of initial temperature ---
cos_perturbation = True
"""
:var: Whether to perturb the temperature initial condition by a function. See :func:`m_equations.Equations.thermal_perturbation`\ .

:vartype: boolean

:meta hide-value:
"""

# --- Thermal perturbation amplitude ---
perturb_ampl 	= 5 # K
"""
:var: Thermal perturbation amplitude :math:`A` in :func:`m_equations.Equations.thermal_perturbation`\ .

:vartype: float, unit K

:meta hide-value:
"""
# --Number of half-cosine waves in lateral direction ---
perturb_freq 	= 1 
"""
:var: Number :math:`N`\ of half-cosine waves in the thermal perturbation formula :func:`m_equations.Equations.thermal_perturbation`\ .

:vartype: float, unit -

:meta hide-value:
"""

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
dt_const = 0.0005 #0.5*kyr*time_scaling

maximum_time_step = True
time_step_scaling = 10 # dt = t / time_step_scaling

scaled_time_step = True
dt_max = 0.1*Myr

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#---------------------------- 6/ RHEOLOGY -----------------------------
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# === PHYSICAL PARAMETERS ===
# --- Gravity acceleratuin ---
g 		= 10.0				# m s^-2
"""
:var: Gravity acceleration.

:vartype: float, unit m s\ :math:`^{-2}`

:meta hide-value:
"""

# --- Density of the domain solid ---
rho_s 	= 3300.0 # kg/m3
"""
:var: Density of the domain solid.

:vartype: float, unit kg m\ :math:`^{-3}`

:meta hide-value:
"""

# --- Density of the liquid below the domain ---
rho_l 	= 1000.0 # kg/m3
"""
:var: Density of the liquid below the domain (subsurface ocean).

:vartype: float, unit kg m\ :math:`^{-3}`

:meta hide-value:
"""

# --- Density of the melt produced by the solid ---
rho_m 	= rho_l # kg/m3
"""
:var: Density of the melt produced by the solid.

:vartype: float, unit kg m\ :math:`^{-3}`

:meta hide-value:
"""

alpha_exp = 1.6e-4	# K^-1

# === PHYSICAL INGREDIENTS ===
# --- Rheology ---
if (str(sys.argv[1]) == "V"):
   elasticity = False

if (str(sys.argv[1]) == "VE"):
   elasticity = True
"""
:var: Whether the material behaves as a Maxwell visco-elastic medium.

:vartype: boolean

:meta hide-value:
"""

plasticity = False
"""
:var: Whether the material yields when the yield stress is reached.

:vartype: boolean

:meta hide-value:
"""

# --- Phase transition at the bottom boundary ----
phase_transition = False
"""
:var: Whether there is a phase transition at the ice-water boundary.

:vartype: boolean

:meta hide-value:
"""

# --- Strength of the phase transition at the ice-water boundary ---
DAL_factor 	= 1e-1 # W/m3
"""
:var: Strength of the phase transition at the ice-water boundary due to the DAL effect

:vartype: float, unit W m\ :math:`^{-3}`

:meta hide-value:
"""

# --- Latent heat of the material ---
Lt 			= 334.0e3 # J/kg
"""
:var: Latent heat of the material.

:vartype: float, unit J kg\ :math:`^{-1}`

:meta hide-value:
"""

# --- Adaptive topography diffusion and topography diffusion factor ---
adaptive_smoothing = False
topo_diff = 1e-8*time_scaling

# --- Rheological parameters ---
# --- VISCOSITY ---
viscosity_type = "composition" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition
"""
:var: Defines which formula for ductile viscosity will be used, see :func:`m_rheology.eta_ductile`\ .

:vartype: string

   * ``"constant"`` for constant viscosity of value :math:`\\eta_0`\ . 
   * ``"temp-dep"`` for temperature-dependent viscosity
   * ``"GK_2001"`` for temperature and stress-dependent viscosity of ice-I
   * ``"composition"`` for composition-dependent viscosity (used for selected benchmarks)

:meta hide-value:
"""

# --- Parameters for constant / temperature-dependent viscosity ---
eta_0 = 1e14	# Pa.s
"""
:var: 
   * Constant viscosity for ``constant`` viscosity type. 
   * Viscosity prefactor for ``temp-dep`` viscosity type. See :func:`m_rheology.eta_ductile`.

:vartype: float, unit Pa s

:meta hide-value:
"""
# --- Activation energy ---
Q_activ = 50e3 	# J/mol
"""
:var: Activation energy for ``temp-dep`` viscosity type.  See :func:`m_rheology.eta_ductile`.

:vartype: float, unit J mol\ :math:`^{-1}`

:meta hide-value:
"""

# --- Grain size ---
d_grain = 1.0e-3
"""
:var: Constant grain size.  See :func:`m_rheology.eta_ductile`.

:vartype: float, unit m

:meta hide-value:
"""

# --- Upper cut-off viscosity ---
eta_max = 1e30
"""
:var: Upper cut-off value for viscoplastic viscosity :math:`\\eta_{vp}`\ .  See :func:`m_rheology.eta_ductile`.

:vartype: float, unit Pa s

:meta hide-value:
"""

# --- Lower cut-off value for plastic viscosity ---
eta_min_plast = 1e-13 #eta_max/1e6
"""
:var: Lower cut-off value for plastic viscosity :math:`\\eta_{p}`\ . See :func:`m_rheology.eta_eff`.

:vartype: float, unit Pa s

:meta hide-value:
"""

# # --- Elasticity ---
# # --- Shear modulus ---
# shear_modulus = 3.52e9
# """
# :var: Shear modulus.

# :vartype: float, unit Pa

# :meta hide-value:
# """

# If elasticity is off, we can compute the new viscosity directly from the strain rate (it will be faster)
# However, the stress formula + VEP iterations would work too
stress_iter_error = 1e-4

# --- Plasticity ---
# --- Angle of internal friction in degrees ---
int_friction_angle = 16.0
int_friction_angle2 = 0.0
"""
:var: Angle of internal friction in degrees.

:vartype: float, unit degrees

:meta hide-value:
"""

# --- Cohesion of an undamaged material ---
cohesion_strong = 400.0
"""
:var: Cohesion of an **undamaged** material. See :func:`m_rheology.cohesion`.

:vartype: float, unit Pa

:meta hide-value:
"""

# --- Cohesion of a fully damaged material ---
cohesion_weak = 20.0
"""
:var: Cohesion of a **fully damaged** material. See :func:`m_rheology.cohesion`.

:vartype: float, unit Pa

:meta hide-value:
"""

yield_stress_max, yield_stress_min = 1e4, 200
"""
:var: Upper and lower cut-off values for the yield stress, respectively.

:vartype: float, float

:meta hide-value:
"""

# --- Value of the plastic strain at which the material starts to accumulate damage ---
eps_strong = 0
"""
:var: Value of the plastic strain at which the material starts to accumulate damage. See :func:`m_rheology.cohesion`.

:vartype: float, unit -

:meta hide-value:
"""

# --- Critical value of the plastic strain beyond which the material is considered as fully damaged ---
eps_weak = 0.1
"""
:var: Critical value of the plastic strain beyond which the material is considered as fully damaged. See :func:`m_rheology.cohesion`.


:vartype: float, unit -

:meta hide-value:
"""

# --- Whether the accumulated plastic strain decreases in time ---
healing = False
"""
:var: Whether the accumulated plastic strain decays in time due to microscopic healing processes. 
      See :func:`m_rheology.plastic_strain_integration`.

:vartype: boolean

:meta hide-value:
"""

# --- Characteristic time scale for the microscopic healing processes ---
healing_timescale = 300*kyr
"""
:var: Characteristic time scale for the microscopic healing processes.
      See :func:`m_rheology.plastic_strain_integration`.

:vartype: float

:meta hide-value:
"""