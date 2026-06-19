3. Thermal convection with phase transition
==============================================


.. figure:: code/demo3/demo3.gif

    Thermal convection in an ice shell with a phase transition at the bottom boundary.
    Note the difference in temperature structure compared to the case with 
    :ref:`closed boundary <1_convection>` . The contours represent the streamlines and their color
    scales the maximum speed at given time step.

-----------------------
Download
-----------------------

Here you can download the :download:`bash script <code/demo3/run_MilleFEuiIle_demo3.sh>`\, the :download:`parameter file <code/demo3/m_parameters_demo3.py>`\ and the :download:`matplotlib script <code/demo3/plot_demo3.py>`\ for the heat flux evolution.
Detailed explanation is given below.

-----------------------------
Parameter file modification
-----------------------------

Compared to the thermal convection in the previous section, we need the bottom boundary to deform,
therefre we prescribe the free surface boundary condition

.. code-block:: python
    :emphasize-lines: 3

    # Boundary conditions for velocity (free_slip, no_slip, free surface, velocity, velocity_x, velocity_y)
    BC_Stokes_problem = [["free_slip"],                # top boundary       (1)
                        ["free_surface"],    # bottom boundary    (2)
                        ["free_slip"],       # left boundary      (3)
                        ["free_slip"]]       # right boundary     (4)

In order to shape the boundary by the phase transition (melting and crystallization),
we turn the phase transition on. The ``DAL_factor`` determines the "strength" of the
*Deguen-Alboussière-Labrosse* effect. The value 10 mW/m :sup:`-3` should result in significant
ocean-shell exchange.

.. code-block:: python
    :emphasize-lines: 2,5

    # --- Phase transition at the bottom boundary ----
    phase_transition = True

    # --- Strength of the phase transition at the ice-water boundary ---
    DAL_factor 	= 1e-2 # W/m3

Finally, instead of temperature-dependent viscosity we prescribe the temperature- and stress-dependent rheology
following Goldsby and Kohlstedt (2001). We choose grain size 1 mm.

.. code-block:: python
    :emphasize-lines: 2,7

    # --- VISCOSITY ---
    viscosity_type = "GK_2001"
    .
    .
    .
    # --- Grain size ---
    d_grain = 1.0e-3

-----------------------------
Visualization of streamlines
-----------------------------

First, in ``plot_demo3.py``, we define points (50 points evenly distributed
along the domain length, placed 2 km above the bottom boundary) whose trajectories will be evaluated.

.. code-block:: python
    :emphasize-lines: 2,3,7

    # --- Starting points for streamlines computation ---
    length = 200e3 # Length of the domain
    num_points = 50
    points = []
    for i in range(1, num_points):
        if (i != int(num_points/2.0)):
            points.append([length/num_points*i, 2e3])


Then, function ``compute_streamlines()`` preforms numerical integration of trajectory of the points
defined above and visualizes the streamlines.