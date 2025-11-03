3. Thermal convection with melting
--------------------------------------

Because of specific formulas for density and constant thermal properties of ice used in the study, the modules with rheology and material properties need to be modified in addition to the parameter file.

.. figure:: gif/convection_melting.gif

    Tidally heated thermal convection in an ice shell following the setup by Tobie et al. (2003). 

----------------------------
Parameter file line by line
----------------------------

The parameter file contains all parameters that can be set. Not all of them 
are needed for every setting, therefore here only the active ones are highlighted.

^^^^^^^^^^^^^^^^^^^^^^^
Output files settings
^^^^^^^^^^^^^^^^^^^^^^^

We start with the directory name 

.. code-block:: python
    :emphasize-lines: 2

    # --- Name of the directory with results ---
    name = "convection"

The directory protection will be turned off, meaning that when lauched repeatedly,
the results will be overwritten each time

.. code-block:: python
    :emphasize-lines: 2

    # Protection from overwriting the directory above
    protect_directory = False

We choose the time units millions of years and output every 25 kyr.

.. There has to be an empty line afterwards
.. code-block:: python
    :emphasize-lines: 3, 6, 9

    # --- In what units the time will be? ---
    # 1.0 - seconds or nondimensional, or yr, kyr or Myr
    time_units = Myr

    # --- String representation of time_units ---
    time_units_string = "Myr"

    # --- Output method for Paraview, HDF5 and tracers ---
    output_frequency = ["time", 25*kyr]

Tracers will be used for partial melt, we need to save them in order to visualize them later.

.. code-block:: python
    :emphasize-lines: 2, 5, 8

    # --- Save tracers into files? ---
    save_tracers = True

    # --- What properties to save? Must be one of the following (order does not matter):
    Tracers_Output = ["melt_fraction"]

    # --- Headers for the columns in the text file for tracers
    Tracers_header = ["Melt fraction (-)"]

To visualize the results in Paraview, we will save temperature, melt fraction, and viscosity.
For initial condition, we will save only temperature.

.. code-block:: python
    :emphasize-lines: 2, 4

    # --- What functions to write into Paraview and HDF5 file ---
    Paraview_Output = ["temperature", "melt_fraction", "viscosity"]

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
    reload_HDF5_step			= 0

    # --- What functions to read from HDF5 file when reloading ---
    reload_HDF5_functions = []

    reload_tracers 	= False

    reload_tracers_step		= 0 # Second column in data_timestamp.dat

    # --- Whether to reset time when reloading ---
    restart_time 	= False

|

^^^^^^^^^^^^^^^^^^^^^^^
Simulation duration
^^^^^^^^^^^^^^^^^^^^^^^

The simulation will be stopped at 20 Myr.

.. code-block:: python
    :emphasize-lines: 2

    # --- Criterion for ending the simulation, e.g. ["time", 1*Myr] or ["step", 1000] ---
    termination_condition = ["time", 20*Myr]

|

^^^^^^^^^^^^^^^^^^^^^^^
Tracers options
^^^^^^^^^^^^^^^^^^^^^^^

Although tracers will be used, the options below are not relevant, because melt is the only mechanism using them, therefore the tracers will be added on the go.

.. code-block:: python

    sm_thickness = 500.0	# m
    om_thickness = 500.0	# m

    integration_method = "RK4" # Euler, RK2, RK4

    tracers_per_cell = 20 

    weight_tracers_by_ratio = False

|

^^^^^^^^^^^^^^^^^^^^^^^
Mesh geometry
^^^^^^^^^^^^^^^^^^^^^^^

The computational domain will be 20 km thick and 40 km long. The mesh will be created
during the inicialization, using crossed geometry of the elements. The basic elements will have 
the longest side of 1 km.

.. code-block:: python
    :emphasize-lines: 2, 5, 8, 11, 14, 17, 20, 23 

    # --- Mesh height ---
    height = 20e3 # m

    # --- Mesh length ---
    length = 40e3 # m

    # --- Whether to read an external mesh ---
    loading_mesh = False

    # --- Mesh to be loaded ---
    mesh_name 	= ""

    # --- Basic mesh resolution if not loading mesh ---
    z_div = 40

    # --- Number of nodes in horizontal direction ---
    x_div = int(z_div*(length/height)) # (keeps aspect ratio 1)

    # --- Method of dividing basic squares into mesh triangle elements. ---
    triangle_types = "crossed" # crossed, left, right, left/right, right/left

    # --- Rectangular mesh refinement ---
    refinement = []

|

^^^^^^^^^^^^^^^^^^^^^^^
Material composition
^^^^^^^^^^^^^^^^^^^^^^^

The ice will consist of pure ice, therefore the list ``materials`` will remain empty.
Other options can be left as they are.

.. code-block:: python
    :emphasize-lines: 1

    materials = []

    empty_cells_allowed = False

    empty_cells_composition = 0

    empty_cells_region = [] 

|

----------------------------
Rheology file modification
----------------------------

-----------------------------------------------
Material properties file modification
-----------------------------------------------

Here we adopt the density of the mixture of ice and melt from Tobie et al. (2003):

.. math::
        
    \rho(T, \chi) = \rho_s(1 - \alpha(T - T_{ref})) + \chi\rho_s(1 - \rho_s / \rho_m)

.. code-block:: python
    :emphasize-lines: 5

    def rho(Temp, composition, xm):
            # .
            # .
            # .
            return rho_s*(1.0 - 1.6e-4*(Temp-T_ref)) + xm*rho_s*(1.0-rho_s/rho_m) 

and constant thermal properties:

.. code-block:: python
    :emphasize-lines: 3,7

    def k(Temp, composition):
            #return 2.3 567.0/Temp
            return 2.3

    def cp(Temp, composition):
            #return 185.0 + 7.037*Temp
            return 2100.0