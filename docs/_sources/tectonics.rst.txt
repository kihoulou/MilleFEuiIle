4. Tectonic extension
----------------------

.. Figure:: gif/3_convection_visc.gif

    Tectonic extension of Europa's ice shell.

----------------------------
Parameter file line by line
----------------------------
.. code-block:: python
    :emphasize-lines: 2

    # --- Name of the directory with results ---
    name = "europa_extension"

The directory protection will be turned off, meaning that when lauched repeatedly,
the results will be overwritten each time

.. code-block:: python
    :emphasize-lines: 2

    # Protection from overwriting the directory above
    protect_directory = False

We choose the time units millions of years and output every 10 kyr.

.. code-block:: python
    :emphasize-lines: 3, 6, 9

    # --- In what units the time will be? ---
    # 1.0 - seconds or nondimensional, or yr, kyr or Myr
    time_units = Myr

    # --- String representation of time_units ---
    time_units_string = "Myr"

    # --- Output method for Paraview, HDF5 and  tracers ---
    output_frequency = ["time", 10*kyr] # e.g. ["steps", 10] or ["time", 100*kyr]

Tracers will be used for elasticity and plasticity, we do not need to save them.
.. code-block:: python
    :emphasize-lines: 2, 5, 8
    # --- Save tracers into files? ---
    save_tracers = False

    # --- What properties to save? Must be one of the following (order does not matter):
    Tracers_Output = []

    # --- Headers for the columns in the text file for tracers
    # --- Up to the user (order corredponding to "stat_output").
    Tracers_header = []

.. code-block:: python
    :emphasize-lines: 2, 4

    # --- What functions to write into Paraview and HDF5 file ---
    Paraview_Output = ["temperature", "viscosity", "plastic_strain", "velocity", "strain_rate_inv", "composition_1"]

    Paraview_Output_Ini = ["temperature"]

To the statistics text file, time (in the units specified above, i.e., Myr)
and the duration of the time step will be saved. We also choose corresponding headers.

.. code-block:: python
    :emphasize-lines: 2, 5
    # --- What values to print in a text file every time step? ---
    stat_output = [ "time", "timestep"]

    # --- Headers for the columns in the text file
    stat_header = ["Time (h)", "dt (s)"]

|

^^^^^^^^^^^^^^^^^^^^^^^
Reloading options
^^^^^^^^^^^^^^^^^^^^^^^
As we define our initial condition and not reload, the following options need to be set.
The number at ``reload_HDF5_step`` and ``reload_tracers_step`` does not matter, but needs
to be defined.


.. code-block:: python
    :emphasize-lines: 2, 5, 8, 11, 13, 15, 18

    # --- Name of the directory from which the data will be reloaded ---
    reload_name 		= ""

    # --- Whether or not to load HDF5 data from data_reload_name/HDF5/data.h5 --- 
    reload_HDF5 		= False

    # ---  Time stamp of the HDF5 file which will be loaded ---
    reload_HDF5_step			= 12

    # --- What functions to read from HDF5 file when reloading ---
    reload_HDF5_functions = ["temperature"]

    reload_tracers 	= False

    reload_tracers_step		= 120 # Second column in data_timestamp.dat

    # --- Whether to reset time when reloading ---
    restart_time 	= False

|

^^^^^^^^^^^^^^^^^^^^^^^
Simulation duration
^^^^^^^^^^^^^^^^^^^^^^^
The simulation will be stopped at 2 Myr.

.. code-block:: python
    :emphasize-lines: 2

    # --- Criterion for ending the simulation, e.g. ["time", 1*Myr] or ["step", 1000] ---
    termination_condition = ["time", 2*Myr]

|

^^^^^^^^^^^^^^^^^^^^^^^
Tracers options
^^^^^^^^^^^^^^^^^^^^^^^

For integration of the tracer trajectories, we choose the 4-th order Runge-Kutta scheme.
Initially, there will be 20 tracers per cell. Weighing by ratio is relevant only in case of multiple
materials, which is now not the case.

.. code-block:: python
    :emphasize-lines: 4, 6

    sm_thickness = 500.0	# m
    om_thickness = 500.0	# m

    integration_method = "RK4" # Euler, RK2, RK4

    tracers_per_cell = 20 

    weight_tracers_by_ratio = False

^^^^^^^^^^^^^^^^^^^^^^^
Mesh geometry
^^^^^^^^^^^^^^^^^^^^^^^

The computational domain will be 25 km thick and 50 km long. The mesh will be created
during the inicialization, using crossed geometry of the elements. The basic elements will have 
the longest side of 1 km.

.. code-block:: python
    :emphasize-lines: 2, 5, 8, 11, 14, 17, 20, 23 

    # --- Mesh height ---
    height = 25e3 # m

    # --- Mesh length ---
    length = 50e3 # m

    # --- Whether to read an external mesh ---
    loading_mesh = False

    # --- Mesh to be loaded ---
    mesh_name 	= "meshes/mesh_25x100km.xml"

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
    refinement = []

|

^^^^^^^^^^^^^^^^^^^^^^^
Material composition
^^^^^^^^^^^^^^^^^^^^^^^

The ice will consist of pure ice, but we will want to track the deformation of the topmost
and bottommost layer of the ice shell. The first layer in the list ``materials`` will be the ice,
the rest the abovementioned layers. The materials will be used solely to track deformation and 
no material or rheological properties will be derived from them.

.. code-block:: python
    :emphasize-lines: 15

    # --- Simple material assignment using "interface", "circle" or "rectangle" geometries ---
    # --- The n-th material overrides the (n-1)th.
    # --- For more complicated geometries, modify "introduce_tracers" function in the m_tracers.py 

    # --- Syntax: [[material_0], [material_1], ...] 
    # Circle:               ["circle",     center_x,   center_y,   radius]
    # Rectangle:            ["rectangle",  left,       right,      bottom,     top]
    # Interface:            ["interface",  "below/above"]

    def interface(x):
        # return 0.02*cos(np.pi*x/length) + 0.2
        return 5e2*cos(2.0*np.pi*x/length) + 3e3

    # --- Leave empty for a single material ---
    materials = [["rectangle", 0, length, 0, height], ["rectangle", 0, length, height- 2e3, height], ["rectangle", 0, length, 0, 2e3]]

    # --- If True, the cells without tracers will be assigned material "default_composition" ---
    # --- Applicable only if the material composition is the only tracer-requiring feature ---
    empty_cells_allowed = False

    empty_cells_composition = 0

    empty_cells_region = [] #  So far rectangle implemented

|

^^^^^^^^^^^^^^^^^^^^^^^
Stokes problem settings
^^^^^^^^^^^^^^^^^^^^^^^

For the discretization of the velocity, we choose Mini elements
The time step will be constant, with the velue defined later in this file. The CFL coefficient will be 
0.5. Since the rheology will be non-Newtonian (stress-dependent viscosity and plasticity), the 
therefore we set the maximym number of Picard iterations and their desired precision.


.. code-block:: python
    :emphasize-lines: 1,3,5,8,11,14

    stokes_elements = "Mini"

    time_step_position = "stokes" #stokes (right after Stokes problem) or "end" (at the end of the time loop)

    time_step_strategy = "constant"

    # --- The CFL parameter ---
    cfl = 0.5

    # ---  Maximum number of Stokes solver Picard iterations ---
    Picard_iter_max 	= 25

    # --- Minimum relative error in velocity field ---
    Picard_iter_error 	= 1e-3

    error_type          = "maximum" # "maximum or integrated"

For the boundary conditions we prescribe free surface at the top boundary,
divergent velocity at the lateral boundaries and influx from the bottom boundary
to compensate the outflux from the sides. The mesh will be moving, and we will
fix the average position of the boundary.

.. code-block:: python
    :emphasize-lines: 2,3,4,5

    # Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
    BC_Stokes_problem = [["free_surface"],                # top boundary       (1)
                        ["velocity", 0.0, 10e3/Myr],    # bottom boundary    (2)
                        ["velocity_x", -10e3/Myr],       # left boundary      (3)
                        ["velocity_x", 10e3/Myr]]       # right boundary     (4)

    mesh_displacement_laplace = "full"  #"full_laplace" or "z_only"

    # Method for correcting velocity field in case of multiple free surface conditions
    stokes_null = "boundary"

^^^^^^^^^^^^^^^^^^^^^^^
Thermal solver settings
^^^^^^^^^^^^^^^^^^^^^^^

To see whether ascending plume penetrating to the stagnant lid changes the temperature at the surface,
we prescribe the radiative boundary condition. At the surface and ice-water boudnary, we prescribe temperature, while the lateral boundaries 
are adiabatic. The reference temperature is set to 273 K. Parameters related to radiative boundary condition
are set to fit Europa.

.. code-block:: python
    :emphasize-lines: 2,5,8,9,10,11,14
        
    # --- Solve heat transfer? ---
    solve_energy_problem = True

    # --- Whether the heat transfer should be solved with nonlinear solver ---
    nonlinear_heat_equation = True

    # --- Boundary condition for heat transfer equation ---
    BC_heat_transfer   = [["radiation", 90],     # top boundary    (1)
                        ["temp", 265.0],     # bottom boundary (2)
                        ["heat_flux", 0.0],  # left boundary   (3)
                        ["heat_flux", 0.0]]  # right boundary  (4)

    # --- Reference temperature ---
    T_ref 	= 265.0 # K

    # --- Insolation parametes ---
    emis 	= 0.97         	# Ice emissivity
    albedo 	= 0.67
    dist_AU = 5.2			# Object's distance from Sun in AU
    insolation = insol_1AU/dist_AU**2

|

In the absence of tidal heating, the melting is not expected to occur. As an initial condition,
a conductive profile will be solved and perturbed by a cosine wave.

.. code-block:: python
    :emphasize-lines: 3,9,19,22,25

    # --- Melting inside the shell ---
    # --- Whether generate partial melt if the temperature of the solid reaches T_melt ---
    internal_melting = False

    # --- Melting temperature of the solid ---
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

    # --- Thermal perturbation amplitude ---
    perturb_ampl 	= 1 # K
    # --Number of half-cosine waves in lateral direction ---
    perturb_freq 	= 1 

|

^^^^^^^^^^^^^^^^^^^^^^^
Time step settings
^^^^^^^^^^^^^^^^^^^^^^^

To mimic Europa's ice shell, we prescribe corresponding gravity and ice physical properties.
Since the density is directly given by a temperature-dependent formula
in ``m_material_properties.py``, the thermal expansivity here is inactive.

.. code-block:: python
    :emphasize-lines: 2,5,8,11

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
    dt_const = 1.0*kyr

    maximum_time_step = False
    time_step_scaling = 25 # dt = t / time_step_scaling

    scaled_time_step = False
    dt_max = 0.1*Myr

|

^^^^^^^^^^^^^^^^^^^^^^^
Rheology
^^^^^^^^^^^^^^^^^^^^^^^

To mimic Europa's ice shell, we prescribe corresponding gravity and ice physical properties.
Since the density is directly given by a temperature-dependent formula
in ``m_material_properties.py``, the thermal expansivity here is inactive.

.. code-block:: python
    :emphasize-lines: 2,5,8,11

    # --- Gravity acceleratuin ---
    g 		= 1.3				# m s^-2

    # --- Density of the domain solid ---
    rho_s 	= 920.0 # kg/m3

    # --- Density of the liquid below the domain ---
    rho_l 	= 1000.0 # kg/m3

    # --- Density of the melt produced by the solid ---
    rho_m 	= rho_l # kg/m3

    alpha_exp = 1.6e-4	# K^-1

    # === PHYSICAL INGREDIENTS ===
    # --- Rheology ---
    elasticity = True

    plasticity = True

    # --- Phase transition at the bottom boundary ----
    phase_transition = True

    # --- Strength of the phase transition at the ice-water boundary ---
    DAL_factor 	= 0.0 # W/m3

    # --- Latent heat of the material ---
    Lt 			= 334.0e3 # J/kg

    # --- Adaptive topography diffusion and topography diffusion factor ---
    adaptive_smoothing = False
    topo_diff = 1e-8*time_scaling

    # --- Rheological parameters ---
    # --- VISCOSITY ---
    viscosity_type = "GK_2001" # constant, temp-dep, GK_2001 (Goldsby and Kohlstedt, 2001) or composition

    # --- Parameters for constant / temperature-dependent viscosity ---
    eta_0 = 1e15	# Pa.s
    # --- Activation energy ---
    Q_activ = 60e3 	# J/mol

    # --- Grain size ---
    d_grain = 1e-3

    # --- Upper cut-off viscosity ---
    eta_max = 1e23

    # --- Lower cut-off value for plastic viscosity ---
    eta_min_plast = eta_max/1e8

    # # --- Elasticity ---
    # If elasticity is off, we can compute the new viscosity directly from the strain rate (it will be faster)
    # However, the stress formula + VEP iterations would work too
    stress_iter_error = 1e-4

    # --- Plasticity ---
    # --- Angle of internal friction in degrees ---
    int_friction_angle = 16.0
    int_friction_angle2 = 0.0

    # --- Cohesion of an undamaged material ---
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
    healing_timescale = 300*kyr
