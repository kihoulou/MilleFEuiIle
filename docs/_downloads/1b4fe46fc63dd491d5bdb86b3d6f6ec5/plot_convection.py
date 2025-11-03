from dolfin import *
import matplotlib.pyplot as plt

import numpy as np
import subprocess
from plot_millefeuille import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size= MPI.size(comm)

name = "europa_melting"

skip_plotting = False

plot_streamlines = False

plot_tracers = False

plot_melt_tracers = True

show_time = False

plot_scalebar = False

plot_axes = True

i_start = 0 

if (skip_plotting == False):
    font_size = 14

    x_label = ["0", "20", "40"]
    x_label_val = [0, 20e3, 40e3]
    x_label_text = "x (km)"

    z_label = ["0", "20"]
    z_label_val = [0, 20e3]
    z_label_text = "z (km)"

    c_label = ["100", "185", "270"]
    c_label_val = [100, 185, 270]
    c_label_text = "Temperature (K)"

    m_label = ["0", "3", "6"]
    m_label_val = [0, 3e-2, 6e-2]
    m_label_text = "Melt fraction (%)"

    # c_label = ["14", "18", "23"]
    # c_label_val = [14, 18, 23]
    # c_label_text = r"log$_{10}$ viscosity (Pa s)"

    x_range = [0, 40e3]
    z_range = [0, 20e3]

    length = 50e3
    points = []
    for i in range(1, 49):
        if (i != 25):
            points.append([length/50.0*i, 15e3])


    labels = [[x_label, x_label_text, x_label_val],
            [z_label, z_label_text, z_label_val],
            [c_label, c_label_text, c_label_val],
            [m_label, m_label_text, m_label_val],]

    # --- Snapshot ---
    # i_HDF5 = 4

    fig = plt.figure(figsize=(7,3))
    ax = plt.subplot2grid(shape=(1, 1), loc=(0, 0))

    # Temp, mesh = initialize(name, i_HDF5)
    # Temp = load(i_HDF5, Temp, name)

    # plot_scalar(ax, fig, Temp, mesh, labels, x_range, z_range, font_size)
    # plt.savefig("img_example.png", dpi = 200, transparent=False, pad_inches=0.3, bbox_inches='tight')

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
            Temp, v, visc, mesh, visc_log, sDG0 = initialize(name, hdf5_step)
            Temp = load(hdf5_step, Temp, name)

            visc_log.assign(project(ln(visc)/ln(10.0), sDG0))

            if (plot_scalebar == True):
                scalebar(ax, 50e3, 0, font_size)

            if (plot_axes == False):
                ax.axis('off')

            if (show_time == True):
                ax.text(0, 25e3 + 2e3, r"$t$ = "+str("{:.1f}".format(ii[2]))+" Myr", fontsize=font_size)

            # --- Plot temperature field ---
            print("\nPlotting snapshot", hdf5_step, "/", int(steps[len(steps) - 1][0]))
            plot_scalar(ax, fig, Temp, mesh, labels, x_range, z_range, font_size)

            # --- Plot tracers ---
            if (plot_tracers == True):
                load_tracers(ax, name, sim_step)

            if (plot_melt_tracers == True):
                load_melt_tracers(ax, name, labels, fig, font_size, sim_step)

            # --- Save figure ---
            if (plot_streamlines == True):
                compute_streamlines(ax, v, points)
                plt.savefig("data_" + name + "/img/temp_str_" + str(hdf5_step).zfill(3) + ".png", dpi = 100, transparent=False, pad_inches=0.3, bbox_inches='tight')
            else:
                plt.savefig("data_" + name + "/img/temp_" + str(hdf5_step).zfill(3) + ".png", dpi = 100, transparent=False, pad_inches=0.3, bbox_inches='tight')
                # plt.savefig("data_" + name + "/img/temp_" + str(hdf5_step).zfill(3) + ".pdf", transparent=False, pad_inches=0.3, bbox_inches='tight')
            plt.close(fig)
            plt.clf()

# --- Make animation using ImageMagick ---
print("\nMaking animation.")
if (plot_streamlines == True):
    subprocess.run(["convert -delay 5 data_" + name + "/img/temp_str_*.png data_" + name + "/anim/3_convection_streamlines.gif"], shell=True) 
else:
    subprocess.run(["convert -delay 5 data_" + name + "/img/temp_*.png data_" + name + "/anim/convection_melting.gif"], shell=True) 

# --- Remove the individual nf figures ---
# print("\nRemoving png figures.")    
# subprocess.run(["rm  data_" + name + "/anim/*.png"], shell=True) 
