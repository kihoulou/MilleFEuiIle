**File tree**
-------------------------------------------------------------------------

When running the code, the following tree is created:

::

    /data_name
    ├── /HDF5         
    │   ├── /meshes
    │   ├── data.h5
    │   ├── data_timestamp.dat
    │   ├── subdomains.pvd, *.pvtu, *.vtu
    │   └── mesh.h5, mesh.xdmf
    ├── /anim
    ├── /img
    ├── /paraview        
    │   ├── /initial_condition
    │   ├── *.h5, *.xdmf
    ├── /source_code
    ├── /tracers
    ├── empty_cells.dat
    ├── picard_iter.dat
    └── statistics.dat

Below, the data in individual directories are described.

^^^^^^^^^^^^^^^^^^^^^^^
/HDF5
^^^^^^^^^^^^^^^^^^^^^^^

* The directory ``meshes`` contains files named ``mesh_*.h5``, which is an HDF5 file containing the mesh for each time output. These files can be used for plotting through *matplotlib* and cannot be opened in Paraview.
* The file ``data.h5`` contains functions defined in ``Paraview_Output`` parameter, but always at least those that are needed for restarting the computation. Data from this file are used for plotting and reloading the simulation.
* The file ``data_timestamp.dat`` contains information about the steps when the data were written (time, simulation step, position in HDF5 file)
* Files starting ``subdomains*`` contain Paraview-readable format of the subdomains, which enables to check whether the boundaries were defined and recognized correctly.
* Files ``mesh.h5`` and ``mesh.xdmf`` contain mesh in paraview-readable format.

|

^^^^^^^^^^^^^^^^^^^^^^^
/anim and /img
^^^^^^^^^^^^^^^^^^^^^^^

These are auxiliary directories that can be used in case of plotting images and making animations from the results.

|

^^^^^^^^^^^^^^^^^^^^^^^
/paraview
^^^^^^^^^^^^^^^^^^^^^^^

* Contains pairs of ``*.xdmf`` and ``*.h5`` for each of the physical quantities specified in ``Paraview_Output`` (again at least those that are needed for restarting the computation).
* The directory  ``initial_condition`` contains quantities that were used as an initial condition and specified in ``Paraview_Output_Ini``

|

^^^^^^^^^^^^^^^^^^^^^^^
/source_code
^^^^^^^^^^^^^^^^^^^^^^^

Contains all parts of the source code needed for reproduction of the results.



^^^^^^^^^^^^^^^^^^^^^^^
/tracers
^^^^^^^^^^^^^^^^^^^^^^^

Contains files with tracer properties for each step of output specified in ``Tracers_Output``.