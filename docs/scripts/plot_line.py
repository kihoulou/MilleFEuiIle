from dolfin import *
import numpy
import math
import os

import matplotlib
matplotlib.use('Agg') #suppress opening

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import ticker

import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.ticker


name = "chaos_uniform_5km_copy"

fig = plt.figure(figsize=(5,2))
size_font = 12

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)

ax0 = plt.subplot2grid(shape=(1, 1), loc=(0, 0)) # row, column, down, right
# ax1 = plt.subplot2grid(shape=(1, 2), loc=(0, 1))

tt = []
cc = []

infile = open("data_" + name +"/cohesion.dat", "r") 
lines = infile.readlines() 

header = True
for line in lines:
    sline = line.split("\t")

    if (header==True):
        header = False
        continue

    tt.append(float(sline[0])/24)
    cc.append(float(sline[1])/1e6)

def plot_data(ax, tt, cc):
    ax.tick_params(axis='y',  colors='black')
    ax.set_ylabel(r"$C$"+" (MPa)", labelpad = 15, fontsize=size_font)
    ax.set_xlabel(r"$t$"+" (days)", labelpad = 15, fontsize=size_font)

    # ax.set_ylim(y_min, y_max)
    # ax.set_xlim(0, 0.3)

   
    # ax.set_xticks(np.linspace(0, 0.3, num= len(xxtics)))
    # ax.set_xticklabels(xxtics,fontsize=size_font)

    # ax.set_yticks(np.linspace(y_min, y_max, num=len(yytics)))
    # ax.set_yticklabels(yytics, fontsize=size_font)

    # ax.tick_params(axis='both', which='major', pad=15)

    ax.plot(tt, cc)

plot_data(ax0, tt, cc)

plt.savefig("data_" + name + "/img_cohesion.png", dpi = 200, transparent=False, pad_inches=0.3, bbox_inches='tight')