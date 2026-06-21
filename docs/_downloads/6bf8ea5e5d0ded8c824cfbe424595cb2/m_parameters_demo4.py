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
name = "demo_convection_comp"
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

# --- Name of the directory from which the data will be reloaded ---
reload_name 		= ""
"""
:var: The directory from which the data will be reloaded when ``reloading_HDF5 = True``.

:vartype: string

:meta hide-value:
"""

# --- Whether or not to load HDF5 data from reload_name/HDF5/data.h5 --- 
reload 		= False
"""
:var: Whether or not to load HDF5 data from ``reload_name/HDF5/data.h5``. 

:vartype: boolean

:meta hide-value:
"""

# ---  Time stamp of the HDF5 file which will be loaded ---
reload_step			= 10
"""
:var: Time stamp of the HDF5 file which will be loaded when ``reload = True``.

:vartype: integer or "last" (reloads the last time step)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------------  TIME SETTINGS -------------------
#----------------------------------------------------------------------

# --- In what units the time will be? ---
# 1.0 - seconds or nondimensional, or yr, kyr or Myr
# --- String representation of time_units ---
time_units_string = "Myr"
"""
:var: String representation of ``time_units``, up to the user. 

:vartype: string, e.g.  ``"-"``, ``"s"``, ``"yr"``, ``"kyr"``, ``"Myr"`` or ``"Gyr"``

:meta hide-value:
"""

# --- Criterion for ending the simulation, e.g., ["time", 1*Myr], ["step", 1000] or ["initial_condition", *]---
termination_condition = ["time", 20*Myr]
""" Time or step criterion for ending the simulation naturally.

:var:  
   * ``["step", 1000]`` exits the time loop after the 100th step
   * ``["time", 100*kyr]`` exits the time loop as soon as 100 kyr is reached 
   * ``["initial_condition", 100*kyr]`` exits the time loop after writing the initial condition


:meta hide-value:
"""

# --- Output method for Paraview, HDF5 and  tracers ---
output_frequency = ["time", 100*kyr] # e.g. ["steps", 10] or ["time", 100*kyr]
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

tracers_output = []
"""
List of strings indicating which properties carried by the tracers will be saved, e.g.::

   tracers_output = ["rank", "dev_stress_xx", "melt_fraction"]

.. table:: Properties to save on tracers
   :widths: auto

   =================  ==============
     String           Description
   =================  ==============
   rank                the process the tracer belongs to
   dev_stress_xx       *xx* component of the deviatoric stress
   dev_stress_xz       *xz* component of the deviatoric stress
   plastic_strain      amount of plastic strain
   original_y          original *z*-coordinate of the tracer
   composition_n       concentration of material *n* (*n* being a number)
   melt_fraction       amount of partial melt on the tracer
   origin              0 if the tracer is original, 1 if added later
   ID                  unique ID of the tracer
   =================  ==============

:meta hide-value:
"""

# --- Headers for the columns in the text file for tracers ---
tracers_header = [] # KEEP EMPTY, will be filled automatically

# --- What functions to write into Paraview and HDF5 file ---
paraview_output = ["temperature", "velocity", "viscosity", "composition_1"]
"""
List of strings indicating which properties will be saved, e.g.::

   paraview_output = ["viscosity", "temperature", "plastic_strain"]

The following are defined: 

======================  ===================================================
String                   Description
======================  ===================================================
``temperature``         :math:`T` (K)
``velocity``            :math:`v` (m s :math:`^{-1}`)
``pressure``            :math:`p` (Pa)
``viscosity``           :math:`\eta` (Pa s)
``viscosity_log``       log\ :math:`_{10}(\eta)` (-)
``eta_v``               :math:`\eta_{ductile}` (Pa s)
``composition_0``       1st material in ``materials[]``
``composition_1``       2nd material in ``materials[]``
``composition_2``       3rd material in ``materials[]``
``tracers``             number of tracers in the cell
``tidal_heating``       :math:`H_{tidal}` (W m :math:`^{-3}`)
``melt_fraction``       :math:`\chi_m` (%)  
``mesh_displacement``   :math:`\\boldsymbol{u}_{mesh}` (m)  
``topography_bottom``   :math:`h_{bot}` (-)  
``topography_top``      :math:`h_{top}` (-)
``iteration_error``     
``plastic_strain``      :math:`\\varepsilon_p` (-)
``yield_stress``        :math:`\sigma_Y` (Pa)
``stress_dev_inv``      :math:`\sigma_{II}` (Pa)
``strain_rate_inv``     :math:`\dot{\\varepsilon}_{II}` (s :math:`^{-1}`)
``yield_function``      1 if yielding, 0 otherwise
``density``             :math:`\\rho` (kg m :math:`^{-3}`)
``z_function``          :math:`Z`, visco-elastic parameter
``shear_modulus``       :math:`G` (Pa)
``cohesion``            :math:`C` (Pa)
======================  ===================================================

.. tip:: A function can be added by passing it to ``__init__`` in ``m_filenames.py``. 
         Follow the procedure for already existing functions.
                    
:meta hide-value:
"""

paraview_output_ini = ["temperature", "composition_1"]
"""
List of strings indicating which properties will be saved at initial condition.

The following are defined: 

======================  ===================================================
Input                   Note
======================  ===================================================
``temperature``         :math:`T` (K)
``topography_top``      :math:`h_{top}` (-)
``topography_bottom``   :math:`h_{bot}` (-)  
``tracers``             number of tracers in the cell
``density``             :math:`\\rho` (kg m :math:`^{-3}`)
``melt_fraction``       :math:`\chi_m` (%)  
======================  ===================================================

:meta hide-value:
"""

# --- What values to print in a text file every time step? Must be one of the following (order does not matter):
# nusselt 	= Nusselt number
# vrms 		= Root mean square velocity
# tracers 	= Number of tracers
# q_top 	= Heat flux over the top boundary
# q_bot 	= Heat flux over the bottom boundary
# time		= Duration of the simulation (hours)
# timestep	= Duration of the time step (seconds)
stat_output = ["q_top", "time", "timestep"]
"""
List of strings indicating which properties will be saved in the text file at each time step.

The following are defined: 

======================  ===================================================
Input                   Note
======================  ===================================================
``nusselt``             :math:`Nu` (-)
``q_top``               :math:`q_{top}` (W m :math:`^{-2}`)
``vrms``                :math:`v_{rms}` (m s :math:`^{-1}`)
``avg_h_bot``           average bottom topography (m)
``h_top_max``           maximum top topography (m)
``time``                wallclock time of the simulation (hours)
``timestep``            one time step duration (seconds)
======================  ===================================================

:meta hide-value:
"""

# --- Headers for the columns in the text file
# --- Up to the user (order corredponding to "stat_output").
stat_header = ["q_top (W/m2)", "Time (h)", "dt (s)"]

monitor_cache = False
"""
:var: If ``True``, returns the number of cache in each time step strored in ``/home/username/.cache/dijitso/lib/``.

:vartype: boolean

.. warning:: Number of caches should stop increasing after a few time steps. If not, FEnICS is making 
             a new cache in each time step, which may lead to crash of the code after some time.

:meta hide-value:
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------ 3/ MESH SETTINGS ----------------------------
#----------------------------------------------------------------------

# --- Mesh height ---
height = 100e3 # m
""" Height of the rectangular mesh.

:vartype: float

:meta hide-value:
"""

# --- Mesh length ---
length = 200e3 # m
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
mesh_name 	= ""
""" Mesh to be loaded.

:var: path to a mesh in **xml** format, e.g. ``mesh_25x100km.xml``

:vartype: string

:meta hide-value:
"""

# --- Basic mesh resolution if not loading mesh ---
z_div = 20
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
refinement = [0, length, 0, height/2, 80e3, 120e3, height/2 + 1e3, 80e3]
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
    """A function that describes the shape of an interface.

    :meta hide-value:.
    """
    return 75e3*exp(-(x - 100e3)**2/1e8)

# --- Leave empty for a single material ---
materials = [["rectangle",  0, length, 0, height], ["interface",  "below"]]
"""
:var: Simple material assignment using "interface", "circle" or "rectangle" geometries.
:vartype: List of lists consisting of the region shape and its dimensions:

   * Circle: ``["circle", center_x, center_y, radius]``
   * Rectangle:  ``["rectangle", left, right, bottom, top]``
   * Interface: ``["interface", "below/above"]``

   For example, in order to prescribe an ice shell with 1-km-thick distinct layer at the surface, enter the following::

      materials = [["rectangle",  0, length, 0, height],  ["rectangle",  0, length, height-1e3, height]]

.. attention:: The n-th material overrides the (n-1)th.
.. tip:: For more complicated geometries, modify ``introduce_tracers`` function in the m_tracers.py 

:meta hide-value:
"""

# --- If True, the cells without tracers will be assigned material "default_composition" ---
# --- Applicable only if the material composition is the only tracer-requiring feature ---
empty_cells_allowed = True
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

empty_cells_region = [["rectangle", 0, length, 80e3, height]] #  So far rectangle implemented
"""
:var: Specifies geometry of the region where the empty cells will be given  ``empty_cells_composition`` even without tracers inside,
      see :func:`m_tracers.Tracers.introduce_tracers`\ .
:vartype: list of lists consisting of the region shape and its dimensions.
   For example ::

      empty_cells_region = [["rectangle",  0, length, 20e3, height]]

.. note:: Only ``"rectangle"`` type is implemented at this point.

:meta hide-value:
"""

initial_topography = False
"""
:var: Whether there will be an initial topography or not.`\ .

:vartype: boolean

:meta hide-value:
"""

h_top_ini = Expression("0.0*sin(120*pi*x[0]/l)", l = length, pi = np.pi, degree = 1)
"""
:var: Expression defining initial top topography. Will be take into account only if ``initial_topography = True``. 

:vartype: Expression

:meta hide-value:
"""

h_bot_ini = Expression("8.5e3*(tanh((x[0]-l/4)/3000) - tanh((x[0] - 3*l/4)/3000))/2.0", l = length, degree = 1)
"""
:var: Expression defining initial bottom topography. Will be take into account only if ``initial_topography = True``. 

:vartype: Expression

:meta hide-value:
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------ 4/ STOKES PROBLEM SETTINGS -------------------
#----------------------------------------------------------------------
stokes_elements = "Mini"

time_step_position = "stokes"# (right after Stokes problem) or "end" (at the end of the time loop)
"""Determines where to compute new time step.

:var: ``"stokes"`` after the Stokes solver or  ``"end"`` at the end of the time loop
:vartype: string

:meta hide-value:
"""

time_step_strategy = "domain"
""" Specifies the method for computing a new time step, see :func:`m_timestep.time_step`\ .

:var:
   
   * ``"domain"`` computes a new time step based on specific criteria in every step of the time loop (see ``m_timestep`` module)
   * ``"convective"`` prescribes a convective time step given by :math:`v_{max}/\\Delta x_{min}` 
   * ``"constant"`` prescribes a constant time step given by the ``dt_const`` parameter

:vartype: string

:meta hide-value:
"""

cfl = 0.5
""" 
:var: The CFL parameter, determining the length of the adaptive time step, see :func:`m_timestep.time_step`\ .
:vartype: float

:meta hide-value:
"""

# ---  Maximum number of Stokes solver Picard iterations (used for heat equation, advection of tracers) ---
Picard_iter_max 	= 10
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

error_type          = "integrated" # "maximum or integrated"
""" Determines the way of calculating the iteration error.

:var: 
   
   * ``"integrated"`` computes the error as :math:`\\, \\int_\\Omega |\\boldsymbol{v} - \\boldsymbol{v_k}|/|\\boldsymbol{v}| \\textrm{d}\\Omega` 
   * ``"maximum"`` computes the error as :math:`\\, \\textrm{max}(\\boldsymbol{v} - \\boldsymbol{v_k}|/| \\boldsymbol{v}|)` 

:vartype: string
   
:meta hide-value:
"""

# Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
BC_Stokes_problem = [["free_slip"],                # top boundary       (1)
                     ["free_slip"],    # bottom boundary    (2)
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

mesh_displacement_laplace = "full"  #"full" or "z_only"
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

integration_method = "RK4" # Euler, RK2, RK4
""" Method for integration of the trajectories of Lagrangian tracers.

:var:  
   * ``"RK2"`` for the 2nd order Runge-Kutta scheme
   * ``"RK4"`` for the 4nd order Runge-Kutta scheme

:vartype: string

:meta hide-value:
"""

tracers_per_cell = 10

weight_tracers_by_ratio = False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#:::::::::::::::::::: HEAT TRANSFER SETTINGS ::::::::::::::::::::::::::
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- Solve heat transfer? ---
solve_energy_problem = True
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

periodic_BCs = False
"""
:var: If ``True``, the values for tempecrature, pressure and velocity at the left and right boundary are identical (i.e., solves the problem on a surface of a cylinder).

:vartype: boolean

.. attention:: Overrides left and right BCs for thermal and Stokes solver.

.. attention:: Tested only without free surface on top/bottom boundary.

:meta hide-value:
"""

# --- Boundary condition for heat transfer equation ---
BC_heat_transfer   = [["temp", 90.0],     # top boundary    (1)
                     ["temp", 265.0],     # bottom boundary (2)
                     ["heat_flux", 0.0],  # left boundary   (3)
                     ["heat_flux", 0.0]]  # right boundary  (4)
"""
:var: Boundary conditions for heat transfer equation in the following order

   * top boundary ``["temp/heat_flux/radiation", value]`` 

   .. note:: 
      If ``"radiation"`` boundary condition is selected, the ``value`` serves as an initial estimate of the top boudnary temperature for the nonlinear solver.

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

# --- Threshold for adding or removing tracers when melting ---
dT_melt_tresh = 4.0
"""
:var: Threshold for adding or removing tracers when melting.

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
heating_model               = "Andrade" # Maxwell, Andrade or none
H_max = 4e-6 # W m^{-3}

# Andrade parameters
alpha_and = 0.2

# --- Find conductive initial condition ---
init_cond_profile = True
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
# --- Gravity acceleration ---
g 		= 1.3				# m s^-2
"""
:var: Gravity acceleration.

:vartype: float, unit m s\ :math:`^{-2}`

:meta hide-value:
"""

# --- Rheology ---
elasticity = False
"""
:var: Whether the material behaves as a Maxwell visco-elastic medium.

:vartype: boolean

:meta hide-value:
"""

G_ice = 3.52e9   # Pa
"""
:var: Shear modulus of ice.

:vartype: float, unit Pa

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
DAL_factor = 1e-3 # W/m3
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
topo_diff = 1e-8

# --- Rheological parameters ---
# --- VISCOSITY ---
viscosity_type = "GK_2001" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition
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
eta_0 = 2e14 # Pa.s
"""
:var: 
   * Constant viscosity for ``constant`` viscosity type. 
   * Viscosity prefactor for ``temp-dep`` viscosity type. See :func:`m_rheology.eta_ductile`.

:vartype: float, unit Pa s

:meta hide-value:
"""

# --- Activation energy ---
Q_activ = 60e3 	# J/mol
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
eta_max = 1e23
"""
:var: Upper cut-off value for viscoplastic viscosity :math:`\\eta_{vp}`\ .  See :func:`m_rheology.eta_ductile`.

:vartype: float, unit Pa s

:meta hide-value:
"""

# --- Lower cut-off value for plastic viscosity ---
eta_min_plast = 1e10 #eta_max/1e6
eta_min = 1e10 #eta_max/1e6
"""
:var: Lower cut-off value for plastic viscosity :math:`\\eta_{p}`\ . See :func:`m_rheology.eta_eff`.

:vartype: float, unit Pa s

:meta hide-value:
"""

# # --- Elasticity ---
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

# --- Cohesion of an undamaged material (Pa) ---
cohesion_strong = 1e6
"""
:var: Cohesion of an **undamaged** material. See :func:`m_rheology.cohesion`.

:vartype: float, unit Pa

:meta hide-value:
"""

# --- Cohesion of a fully damaged material ---
cohesion_weak = 0.0
"""
:var: Cohesion of a **fully damaged** material. See :func:`m_rheology.cohesion`.

:vartype: float, unit Pa

:meta hide-value:
"""

yield_stress_max, yield_stress_min = 1e8, 0.1
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
eps_weak = 0.2
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
healing_timescale = 100*kyr
"""
:var: Characteristic time scale for the microscopic healing processes.
      See :func:`m_rheology.plastic_strain_integration`.

:vartype: float

:meta hide-value:
"""
