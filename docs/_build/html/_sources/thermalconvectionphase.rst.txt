2. Thermal convection with phase transition
==============================================


.. figure:: gif/2a_convection.gif

    Thermal convection in an ice shell with a phase transition at the bottom boundary.
    Note the difference in temperature structure compared to the case with 
    :ref:`closed boundary <1_convection>` .


-----------------------
Parameter file
-----------------------

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

Function ``compute_streamlines()`` allows to visualize streamlines

.. figure:: gif/img_str_088.png

    Thermal convection with a phase transition. Open streamlines at the ice-water interface
    demonstrate direct material exchange between the subsurface ocean and the ice shell.

-----------------------
Download
-----------------------

Here you can download the complete :download:`parameter file <code/MilleFEuiIle.zip>`\  and the :download:`plotting script <scripts/plot_convection.py>`\ for the animation.
