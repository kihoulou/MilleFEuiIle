from dolfin import *
import matplotlib.pyplot as plt

import numpy as np
import subprocess
from plot_MilleFEuiIle import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

name = "demo_convection_a1"
"""
:var: Name of the directory from which the results will be read.


:vartype: string

:meta hide-value:
"""

skip_plotting = False
"""
:var: Whether to plot the data ``True`` or to skip plotting and directly make an animation from already existing figures ``False``


:vartype: boolean

.. warning:: If the animation memory requirements are too high for ``convert -delay``, making of the gif can fail. In that case, some time steps can be skipped,
             the resolution can be lowered, or some external tools may be used to make the animation.


:meta hide-value:
"""

plot_HDF5_data = True
"""
:var: Whether to plot the HDF5 data from the simulation - temperature, velocity, viscosity, etc. For example, if ``False``, enables the tracers to be plotted alone.


:vartype: string

:meta hide-value:
"""

plot_streamlines = False
"""
:var: Whether to plot streamlines based on the velocity field.


:vartype: boolean

:meta hide-value:
"""

plot_comp_tracers = False
"""
:var: Whether to plot streamlines based on the velocity field.


:vartype: boolean

:meta hide-value:
"""

# --- IDs of the tracers to be plotted ---
plot_trajectories = [["a", "white"], ["b", "aqua"], ["c", "violet"], ["d", "yellow"], ["e", "lawngreen"]]
"""
:var: Whether to plot streamlines based on the velocity field.

:vartype: boolean


.. hint::
    For example, to plot the tracers with IDs from "a" to "e" (defined in the parameter file):

    .. code-block:: python
        :emphasize-lines: 1

        plot_trajectories = [["a", "white"], ["b", "aqua"], ["c", "violet"], ["d", "yellow"], ["e", "lawngreen"]]


:meta hide-value:
"""

tail = "comet"
"""
:var: What type of tail to use for plotting the trajectory of the tracers.


:vartype: ``"full"``, ``"comet"`` or ``"none"``

|
    
The option ``"full"``  shows the full trajectory, but the image soon becomes unclear.

.. figure:: code/other/demo_full.gif

    
The option ``"comet"``  creates a comet-like tail which increases the clarity of the figure,
as the tail length is given by the ``tail_length`` parameter (here set to 500  steps).

.. figure:: code/other/demo_comet.gif


:meta hide-value:
"""

tail_length = 500 # steps
"""
:var: Length of the tail behind the tracer in terms of time steps.


:vartype: integer

:meta hide-value:
"""

plot_melt_tracers = False
"""
:var: Whether to plot melt tracers.

:vartype: boolean

|

Simulation of tidally-induced internal melting (see Demo 5) with melt content shown on tracers.

.. figure:: code/demo5/demo5.gif


:meta hide-value:
"""

show_time = True
"""
:var: Whether to show time in the simulation. Units and position in the figure can be adjusted on the following line

.. code-block:: python
    :emphasize-lines: 2

    if (show_time == True):
        ax.text(0, 1.1*z_label_val[-1], r"$t$ = "+str("{:.1f}".format(ii[2]))+" Myr", fontsize = font_size)

:vartype: boolean

:meta hide-value:
"""

plot_scalebar = False
"""
:var: Plots a rectangular scale bar inside the domain.

:vartype: boolean

.. figure:: code/other/scalebar.png


:meta hide-value:
"""

plot_axes = True
"""
:var: Puts labels on the *x*- and *y*-axis.

:vartype: boolean

.. figure:: code/other/axes.png

:meta hide-value:
"""

i_start = 0 
"""
:var: From which step to plot? (So that the script does not need to replot already plotted data.)

:vartype: Integer

:meta hide-value:
"""

def define_labels():
    """
    :var: Defines labels for x and y axis, and for the scale bars.

    :vartype: List of lists

    .. hint::

        Firs list defines what will be writen on the axis 

        .. code-block:: python
            :emphasize-lines: 1
        
            x_label = ["0", "100", "200"]

        Second line defines what will be the corresponding numerical value 

        .. code-block:: python
            :emphasize-lines: 1

            x_label_val = [0, 100e3, 200e3]

        Third line defines what will be the label of the axis.

        .. code-block:: python
            :emphasize-lines: 1

            x_label_text = "x (km)"


    :meta hide-value:
    """

    x_label = ["0", "100", "200"]
    x_label_val = [0, 100e3, 200e3]
    x_label_text = "x (km)"

    z_label = ["0", "100"]
    z_label_val = [0, 100e3]
    z_label_text = "z (km)"

    # --- Values for temperature color bar ---
    c_label = ["90", "185", "265"]
    c_label_val = [90, 185, 265]
    c_label_text = "Temperature (K)"

    # --- Values for viscosity color bar ---
    # c_label = ["14", "18", "23"]
    # c_label_val = [14, 18, 23]
    # c_label_text = r"log$_{10}$ viscosity (Pa s)"

    # --- Values for melt fraction color bar ---
    m_label = ["0", "3", "6"]
    m_label_val = [0, 3e-2, 6e-2]
    m_label_text = "Melt fraction (%)"

    labels = [[x_label, x_label_text, x_label_val],
            [z_label, z_label_text, z_label_val],
            [c_label, c_label_text, c_label_val],
            [m_label, m_label_text, m_label_val]]

    return labels

def define_streamlines():
    """
    :var: Defines the coordinates of the points from which the streamlines will be calculated.

    :vartype: List of lists 

    :meta hide-value:
    """

    # --- Starting points for streamlines computation ---
    length = 50e3
    points = []
    for i in range(1, 49):
        if (i != 25):
            points.append([length/50.0*i, 15e3])

    return points

font_size = 14

# --- Define labels ---
labels = define_labels()

# --- Define points for streamline computation ---
points = define_streamlines()

x_range = [labels[0][2][0], labels[0][2][-1]]
# z_range = [labels[1][2][0], labels[1][2][-1]]
z_range = [-1e3, labels[1][2][-1]]

if (skip_plotting == False):
    def run_code():
        """
        Runs the plotting loop.

        :meta hide-value:
        """
        fig = plt.figure(figsize=(7,3))
        ax = plt.subplot2grid(shape=(1, 1), loc=(0, 0))

        # --- Animation ---
        infile = open("data_" + name + "/HDF5/data_timestamp.dat", "r") 
        lines = infile.readlines() 

        every_n = 1
        steps = []

        header = True
        for line in lines:
            sline = line.split("\t\t")

            if (header==True):
                header = False
                continue

            steps.append([int(sline[0]), int(sline[1]), float(sline[3])])

        steps = steps[0::every_n]

        for ii in steps:

            if (ii[0] > i_start):
                hdf5_step = ii[0]
                sim_step = ii[1]

                fig = plt.figure(figsize=(7,3))
                ax = plt.subplot2grid(shape=(1, 1), loc=(0, 0))

                ax.spines[['top']].set_visible(False)

                # ---- In each time step, initialize mesh, funciton spaces and functions ---
                # (necessary when the mesh changes shape)
                mesh, sCG2, vCG1, sDG0 = initialize(name, hdf5_step)

                Temp = Function(sCG2)
                v = Function(vCG1)

                Temp, v  = load(name, hdf5_step,
                                "temperature", Temp,
                                "velocity", v
                                )

                if (plot_scalebar == True):
                    scalebar(ax, 190e3, 2e3, font_size)

                if (plot_axes == False):
                    ax.axis('off')

                if (show_time == True):
                    ax.text(0, 1.1*labels[1][2][-1], r"$t$ = "+str("{:.1f}".format(ii[2]))+" Myr", fontsize=font_size)

                

                # --- Plot tracers ---
                if (len(plot_trajectories) > 0):
                    for j in range(0, len(plot_trajectories)):
                        plot_trajectory(ax, name, plot_trajectories[j], sim_step, tail, tail_length)

                # --- Plot temperature field ---
                print("\nPlotting snapshot", hdf5_step, "/", int(steps[len(steps) - 1][0]))
                plot_scalar(ax, fig, Temp, mesh, labels, x_range, z_range, font_size, "smooth", plot_HDF5_data)

                if (plot_comp_tracers == True):
                    load_comp_tracers(ax, name, sim_step)

                if (plot_melt_tracers == True):
                    load_melt_tracers(ax, name, labels, fig, font_size, sim_step)

                # --- Save figure ---
                if (plot_streamlines == True):
                    compute_streamlines(ax, v, points)
                    plt.savefig("data_" + name + "/img/temp_str_" + str(hdf5_step).zfill(3) + ".png",
                                dpi = 250, transparent = False, pad_inches = 0.3, bbox_inches = 'tight')
                else:
                    plt.savefig("data_" + name + "/img/temp_" + str(hdf5_step).zfill(3) + ".png",
                                dpi = 250, transparent = False, pad_inches = 0.3, bbox_inches = 'tight')
                plt.close(fig)
                plt.clf()

        # --- Make animation using ImageMagick ---
        print("\nMaking animation.")
        subprocess.run(["convert -delay 10 data_" + name + "/img/temp_*.png data_" + name + "/anim/demo_full.gif"], shell=True) 

if __name__ == "__main__":
    run_code()