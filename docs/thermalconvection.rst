1. Thermal convection 
----------------------

.. figure:: gif/convection.gif

    Thermal convection in an ice shell.

.. rubric:: Parameter file

We choose the time units millions of years and output every 10 kyr.

.. There has to be an empty line afterwards
.. code-block:: python
    :emphasize-lines: 3, 6, 9

    # --- In what units the time will be? ---
    # 1.0 - seconds or nondimensional, or yr, kyr or Myr
    time_units = Myr

    # --- String representation of time_units ---
    time_units_string = "Myr"

    # --- Output method for Paraview, HDF5 and tracers ---
    output_frequency = ["time", 10*kyr] # e.g. ["steps", 10] or ["time", 100*kyr]

.. rubric:: Download 

Here you can download the complete :download:`parameter file <code/MilleFEuiIle.zip>`\  and the :download:`plotting script <scripts/plot_convection.py>`\ for the animation.
