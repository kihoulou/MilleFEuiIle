.. MilleFEuiIle documentation master file, created by
   sphinx-quickstart on Sat Mar 15 14:02:03 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. MilleFEuiIle  documentation
.. =============================

.. .. image:: _static/mf.png

.. Welcome to the *MilleFEuiIle*'s documentation!

The scope of *MilleFEuiIle*
----------------------------

*MilleFEuiIle* [mil·fœj], a short for **M**\ ulti-\ **I**\ nstrumenta\ **L** **L**\ agrangian-\ **E**\ ulerian **FE**\ niCS-\ **U**\ sing
**I**\ ce-\ **I** **L**\ ayer **E**\ xplorer, is a computational code built 
upon the open-source finite-element python library `FEniCS <https://fenicsproject.org/>`_, developed to solve thermo-mechanical 
evolution of ice shells of icy bodies. Combining finite element method and Lagrangian markers, 
*MilleFEuiIle* can be used to address a variety of problems in 2D Cartesian geometry, such as


* thermal and thermo-compositional convection
* evolution of surface and ice-ocean interface
* phase transition at the ice-ocean interface
* tectonic deformation with visco-elasto-plastic rheology
* melt generation and transport inside the ice shell (using Rayleigh-Taylor approximation)

|

Installation and running
-------------------------

You can download the most recent version of MilleFEuiIle :download:`here <code/MilleFEuiIle.zip>`\ . It 
is built on `FEniCS 2019.2.0 <https://fenicsproject.org/download/archive/>`_, which is currently installed at the `cluster <https://geo.mff.cuni.cz/pocitace.htm>`_ of the Department of Geophysics.

To launch at *N* cores, type into the command line

.. code-block:: bash
   :emphasize-lines: 1

   mpirun -n N main.py > log.out 2> log_e.out&

You can also use the predefined **bash script**, where the number of cores and the name of the current parameter file can be preset (e.g., ``m_parameters_europa.py``).
The script renames the main file and the terminal output accordingly (``main_europa.py``, ``europa.out`` and ``europa_e.out``),
and runs the simulation on background

.. code-block:: bash
   :emphasize-lines: 1

   bash run_code &

|

Visualization
-------------------------
For visualization of the mesh data, use `Paraview <https://www.paraview.org/>`_, or the example matplotlib scripts in :ref:`demos`.



|

Citation
-------------------------

Thermal convection, free surface, **phase boundary evolution**:

* **Kihoulou, M.**, Čadek, O., Kvorka, J., Kalousová, K., Choblet, G. and Tobie, G. (2023). Topographic response to ocean heat flux anomaly on the icy moons of Jupiter and Saturn. *Icarus* **391**, 115337, `https://doi.org/10.1016/j.icarus.2022.115337 <https://doi.org/10.1016/j.icarus.2022.115337>`_
 

Thermal convection, free surface, **tectonic deformation**:

* **Kihoulou, M.**, Choblet, G., Tobie, G., Kalousová, K. and Čadek, O. (2025). Subduction-like process in Europa’s ice shell triggered by enhanced eccentricity periods. *Science Advances* **11**, eadq8719, `https://www.science.org/doi/10.1126/sciadv.adq8719 <https://www.science.org/doi/10.1126/sciadv.adq8719>`_

|

Support
--------------------------

Information on the methodology and the benchmarks can be found in my `dissertation thesis. <https://geo.mff.cuni.cz/theses/2025-Kihoulou-PhD.pdf>`_

This page was last updated on |today| and is still under development. With any questions, issues or ideas, please contact 
me at `martin.kihoulou@mff.cuni.cz <mailto:martin.kihoulou@mff.cuni.cz>`_\ .

|

Explore
--------------------------
.. toctree::
   :maxdepth: 3
   :includehidden:

   demos
   modules