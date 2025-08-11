from dolfin import *
import matplotlib.pyplot as plt

import numpy as np
import subprocess
from plot_millefeuille import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size= MPI.size(comm)

name = "convection_comp_1003"

skip_plotting = False

plot_streamlines = False

plot_tracers = True

if (skip_plotting == False):
    font_size = 10

    x_label = ["0", "100", "200"]
    x_label_val = [0, 100e3, 200e3]
    x_label_text = "x (km)"

    z_label = ["0", "50", "100"]
    z_label_val = [0, 50e3, 100e3]
    z_label_text = "z (km)"

    c_label = ["90", "180", "265"]
    c_label_val = [90, 180, 265]
    c_label_text = "Temperature (K)"

    x_range = [0, 200e3]
    z_range = [0, 100e3]

    length = 200e3
    points = []
    for i in range(1, 49):
        if (i != 25):
            points.append([length/50.0*i, 15e3])


    labels = [[x_label, x_label_text, x_label_val],
            [z_label, z_label_text, z_label_val],
            [c_label, c_label_text, c_label_val]]

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
    time = []

    header = True
    for line in lines:
        sline = line.split("\t\t")

        if (header==True):
            header = False
            continue

        steps.append([int(sline[0]), int(sline[1])])
        time.append(float(sline[2]))

    steps = steps[0::every_n]

    for ii in steps:
        hdf5_step = ii[0]
        sim_step = ii[1]

        fig = plt.figure(figsize=(7,3))
        ax = plt.subplot2grid(shape=(1, 1), loc=(0, 0))

        # ---- In each time step, initialize mesh, funciton spaces and functions ---
        # (necessary when the mesh changes shape)
        Temp, v, mesh = initialize(name, hdf5_step)
        Temp, v = load(hdf5_step, Temp, v, name)

        # --- Plot temperature field ---
        print("\nPlotting snapshot", hdf5_step, "/", int(steps[len(steps) - 1][0]))
        plot_scalar(ax, fig, Temp, mesh, labels, x_range, z_range, font_size)

        # --- Plot tracers ---
        if (plot_tracers == True):
            load_tracers(ax, name, sim_step)

        # --- Save figure ---
        if (plot_streamlines == True):
            compute_streamlines(ax, v, points)
            plt.savefig("data_" + name + "/img/img_str_" + str(hdf5_step).zfill(3) + ".png", dpi = 100, transparent=False, pad_inches=0.3, bbox_inches='tight')
        else:
            plt.savefig("data_" + name + "/img/img_" + str(hdf5_step).zfill(3) + ".png", dpi = 100, transparent=False, pad_inches=0.3, bbox_inches='tight')
        plt.close(fig)
        plt.clf()

# --- Make animation using ImageMagick ---
print("\nMaking animation.")
if (plot_streamlines == True):
    subprocess.run(["convert -delay 5 data_" + name + "/img/img_str_*.png data_" + name + "/anim/3_convection_streamlines.gif"], shell=True) 
else:
    subprocess.run(["convert -delay 5 data_" + name + "/img/img_*.png data_" + name + "/anim/3_convection.gif"], shell=True) 

# --- Remove the individual nf figures ---
# print("\nRemoving png figures.")    
# subprocess.run(["rm  data_" + name + "/anim/*.png"], shell=True) 