7. Reloading a simulation
-------------------------------------------------------------------------

During reloading, the following data are read

* topography of the top and bottom boundary (if the respective free surface bounday conditions are prescribed)
* mesh-carried quantities, e.g., temperature
* tracer-carried quantities, e.g., melt fraction, plastic strain, etc.

^^^^^^^^^^^^^^^^^^^^^^^
How to reload
^^^^^^^^^^^^^^^^^^^^^^^

First, we need to turn on the reloading

.. code-block:: python
    :emphasize-lines: 1

    reload = True

then specify from what directory 

.. code-block:: python
    :emphasize-lines: 1

    reload_name = "directory"

and which time step 

.. code-block:: python
    :emphasize-lines: 1

    reload_step	= 10

which refers to the ``HDF5_step`` column in the ``data_directory/HDF5/data_timestamp.dat`` file. The time of the first time step can be 
either read from the same file, or restarted to 0, which can be chosen through 

.. code-block:: python
    :emphasize-lines: 1

    restart_time = False

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Reading reloaded data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. If needed, temperature is read from ``data_directory/HDF5/data.h5``.

2. If needed, top and bottom topographies are read from ``data_directory/HDF5/data.h5``.
    Mesh is deformed to accomodate the top and bottom topography. The mesh at this point should
    correspond to the one found in ``data_directory/HDF5/mesh/mesh_*.h5``, however, being rectangular in the beginning,
    it was more straightforward to define the boundaries.
3. If needed, the tracers are loaded from ``data_directory/tracers/step_*.dat``.
    Loading of tracers can take several minutes, since each tracer (often hundreds of thousands) has to be assigned to a specific cell (typically hundreds to thousands). 
    Progress is returned every ten percent.

    After being loaded, FEniCS functions for composition, melt fraction, plastic strain and elastic stress are read, if needed.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
