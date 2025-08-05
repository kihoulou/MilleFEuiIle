from dolfin import *
import matplotlib.pyplot as plt

import numpy as np
import subprocess
from plot_millefeuille import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size= MPI.size(comm)

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

name = "convection"

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

    steps.append(int(sline[0]))
    time.append(float(sline[2]))

steps = steps[0::every_n]

for i in steps:
    fig = plt.figure(figsize=(7,3))
    ax = plt.subplot2grid(shape=(1, 1), loc=(0, 0))

    # ---- In each time step, initialize mesh, funciton spaces and functions ---
    # (necessary when the mesh changes shape)
    Temp, v, mesh = initialize(name, int(i))
    Temp, v = load(i, Temp, v, name)

    # --- Plot temperature field ---
    print("\nPlotting snapshot", int(i), "/", int(steps[len(steps) - 1]))
    plot_scalar(ax, fig, Temp, mesh, labels, x_range, z_range, font_size)

    # --- Plot streamlines ---
    # plot_streamlines(ax, v, points)

    # --- Save figure ---
    plt.savefig("data_" + name + "/anim/img_" + str(int(i)).zfill(3) + ".png", dpi = 100, transparent=False, pad_inches=0.3, bbox_inches='tight')
    plt.close(fig)
    plt.clf()

# --- Make animation using ImageMagick ---
print("\nMaking animation.")    
subprocess.run(["convert -delay 5 data_" + name + "/anim/img_*.png data_" + name + "/anim/1_convection.gif"], shell=True) 

# --- Remove the individual nf figures ---
print("\nRemoving png figures.")    
subprocess.run(["rm  data_" + name + "/anim/*.png"], shell=True) 