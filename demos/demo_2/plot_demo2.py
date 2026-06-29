
# --- Python modules ---
import matplotlib.pyplot as plt 

# --- MilleFEuiIle modules ---
from plot_MilleFEuiIle import *

fig, ax0 = plt.subplots(1, 1, layout='constrained', figsize=(6, 4)) 

# --- Define the directories ---
directories = [
    "data_demo2_convection_1e13Pas/statistics.dat",
    "data_demo2_convection_1e14Pas/statistics.dat",
    "data_demo2_convection_1e15Pas/statistics.dat",
    "data_demo2_convection_1e16Pas/statistics.dat"]

# --- Define the labels ---
labelnames = [r"$10^{13}$ Pa s",
              r"$10^{14}$ Pa s",
              r"$10^{15}$ Pa s", 
              r"$10^{16}$ Pa s"]

# --- Define the line colors ---
colors = ["red", "orange", "forestgreen", "navy"]

# --- Define the line styles ---
styles = ["solid", "solid", "solid", "dashed"]

# --- Set the x-axis limit ---
ax0.set_xlim(0, 10)

# --- Plot the third column in the file ---
column = 2
plot_time(ax0, directories, column, labelnames, colors, styles)

# --- Set the x- and y- axes labels ---
ax0.set_ylabel(r"Heat flux (mW m$^{-2}$)", labelpad = 5)
ax0.set_xlabel("Time (Myr)", labelpad = 5)

# --- Plot the legend ---
ax0.legend(loc="upper left", ncol = 1)

# --- Save the figure ---
plt.savefig("demo2.png", bbox_inches = "tight", dpi = 250)
