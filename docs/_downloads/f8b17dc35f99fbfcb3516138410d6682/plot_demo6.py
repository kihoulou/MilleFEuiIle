# --- Python modules ---
from dolfin import *
import matplotlib.pyplot as plt
import subprocess
import numpy as np
import sys

# --- MilleFEuiIle modules ---
from plot_MilleFEuiIle import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

name = "demo6_extension2"

skip_plotting = False

plot_HDF5_data = True

plot_streamlines = False

plot_comp_tracers = [[2, "orange", 0.25], [1, "lightskyblue", 0.25]]

plot_trajectories = []
tail = "comet"
tail_length = 500 # steps

plot_melt_tracers = False

show_time = True

plot_scalebar = True

plot_axes = False

try:
    i_start = int(sys.argv[1])
except:
    i_start = 0

if (skip_plotting == False):
    font_size = 14

    x_label = ["0", "25", "50"]
    x_label_val = [0, 25e3, 50e3]
    x_label_text = "x (km)"

    z_label = ["0", "25"]
    z_label_val = [0, 25e3]
    z_label_text = "z (km)"

    # --- Values for temperature color bar ---
    # c_label = ["90", "185", "265"]
    # c_label_val = [90, 185, 265]
    # c_label_text = "Temperature (K)"

    # --- Values for viscosity color bar ---
    c_label = ["14", "18", "23"]
    c_label_val = [14, 18, 23]
    c_label_text = r"log$_{10}$ viscosity (Pa s)"

    # --- Values for melt fraction color bar ---
    m_label = ["0", "3", "6"]
    m_label_val = [0, 3e-2, 6e-2]
    m_label_text = "Melt fraction (%)"

    # --- Starting points for streamlines computation ---
    length = 200e3
    num_points = 50
    points = []
    for i in range(1, num_points):
        if (i != int(num_points/2.0)):
            points.append([length/num_points*i, 2e3])

    x_range = [x_label_val[0], x_label_val[-1]]
    z_range = [z_label_val[0], 27e3]

    labels = [[x_label, x_label_text, x_label_val],
            [z_label, z_label_text, z_label_val],
            [c_label, c_label_text, c_label_val],
            [m_label, m_label_text, m_label_val],]

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

        if (ii[0] >= i_start):
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
            visc = Function(sDG0)
            visc_log = Function(sDG0)

            Temp, v, visc  = load(name, hdf5_step,
                            "temperature", Temp,
                            "velocity", v,
                            "viscosity", visc
                            )

            visc_log.assign(project(ln(visc)/ln(10), sDG0))

            if (plot_scalebar == True):
                scalebar(ax, 50e3, 10e3, font_size)

            if (plot_axes == False):
                ax.axis('off')

            if (show_time == True):
                ax.text(0, 1.1*z_label_val[-1], r"$t$ = "+str("{:.1f}".format(ii[2]))+" Myr", fontsize=font_size)

             # --- Plot tracer trajectories ---
            if (len(plot_trajectories) > 0):
                for j in range(0, len(plot_trajectories)):
                    plot_trajectory(ax, name, plot_trajectories[j], sim_step, tail, tail_length)

            # --- Plot temperature field ---
            print("\nPlotting snapshot", hdf5_step, "/", int(steps[len(steps) - 1][0]))
            plot_scalar(ax, fig, visc_log, mesh, labels, x_range, z_range, font_size, "flat", plot_HDF5_data)

            # --- Plot composition tracers ---
            if (len(plot_comp_tracers) > 0):
                load_comp_tracers(ax, name, hdf5_step, plot_comp_tracers)

            ax.text(27e3, -2e3, "Oceanic material", color = 'mediumblue', fontsize=font_size-2)
            ax.text(27e3, 27e3, "Near-surface material", color = 'darkorange', fontsize=font_size-2)

            # --- Plot melt tracers ---
            if (plot_melt_tracers == True):
                load_melt_tracers(ax, name, labels, fig, font_size, sim_step)

            # --- Save figure ---
            if (plot_streamlines == True):
                vel_max = MPI.max(mesh.mpi_comm(), np.abs(v.vector().get_local()).max()) 
                compute_streamlines(ax, v, points, vel_max, fig, font_size, i)

            plt.savefig("data_" + name + "/img/visc_" + str(hdf5_step).zfill(3) + ".png",
                            dpi = 100, transparent = False, pad_inches = 0.3, bbox_inches = 'tight')
            plt.close(fig)
            plt.clf()

# --- Make animation using ImageMagick ---
print("\nMaking animation.")
subprocess.run(["convert -delay 5 data_" + name + "/img/visc_*.png data_" + name + "/anim/demo6.gif"], shell=True) 