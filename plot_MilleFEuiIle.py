from dolfin import *
import matplotlib.tri as tri
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.colors as colors

comm = MPI.comm_world
rank = MPI.rank(comm)
size= MPI.size(comm)

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def initialize(name, i):
    mesh = Mesh()

    mesh_file = HDF5File(comm, "data_" + name + "/HDF5/meshes/mesh_" + str(i) +".h5", "r")
    mesh_file.read(mesh, "/mesh", True)
    mesh_file.close()

    sDG0_elem = FiniteElement("DG", mesh.ufl_cell(), 0)
    sDG0 = FunctionSpace(mesh, sDG0_elem)

    sCG2_elem = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
    sCG2 = FunctionSpace(mesh, sCG2_elem)

    vCG1_elem = VectorElement("Lagrange", mesh.ufl_cell(), 1)
    vCG1 = FunctionSpace(mesh, vCG1_elem)

    return mesh, sCG2, vCG1, sDG0

def load(name, i, *args):
    l = []
    files_hdf5_in={}
    files_hdf5_in['data'] = HDF5File(comm, "data_" + name + "/HDF5/data.h5", "r")

    for ii in range(0, len(args), 2):
        dataset = "/" + args[ii] + "/vector_%d"% (i)
        files_hdf5_in['data'].read(args[ii + 1], dataset)
        l.append(args[ii + 1])

    return l

def plot_scalar(ax, fig, function, mesh, labels, x_range, z_range, font_size, type, plot_HDF5_data):
    # Axes range
    ax.set_ylim(z_range[0],z_range[1])
    ax.set_xlim(x_range[0],x_range[1])

    x_label, x_label_text, x_label_val = labels[0][0], labels[0][1], labels[0][2]
    z_label, z_label_text, z_label_val = labels[1][0], labels[1][1], labels[1][2]
    c_label, c_label_text, c_label_val = labels[2][0], labels[2][1], labels[2][2]
    
    ax.set_aspect("equal")
    mesh_function = function.function_space().mesh()
    C = function.compute_vertex_values(mesh)

    if (plot_HDF5_data == True):
        if (type == "flat"):
            values = []
            for j in range(mesh.num_cells()):
                centroid = Cell(mesh, j).midpoint()
                values.append(function(Point(centroid.x(), centroid.y())))
            
            tpc = ax.tripcolor(mesh2triang(mesh_function), values, shading='flat', vmin = c_label_val[0], vmax = c_label_val[-1], cmap = "jet_r")

        if (type == "smooth"):
            tpc = ax.tripcolor(mesh2triang(mesh_function), C, shading='gouraud',  vmin = c_label_val[0], vmax = c_label_val[-1], cmap = "RdYlBu_r")

    # x-label settings
    ax.set_xlabel(x_label_text, fontsize = font_size)
    ax.set_xticks(x_label_val)
    ax.set_xticklabels(x_label, fontsize = font_size)

    # y-label settings
    ax.set_ylabel(z_label_text, fontsize = font_size)
    ax.set_yticks(z_label_val)
    ax.set_yticklabels(z_label, fontsize = font_size)

    ax.tick_params(pad=10) 

    # --- Colorbar ---
    if (plot_HDF5_data == True):
        fig.subplots_adjust(right=0.9)
        cax = ax.inset_axes([1.05, 0, 0.03, 1])
        cbar = fig.colorbar(tpc, cax=cax, orientation='vertical', ticks = c_label_val)

        cbar.ax.set_yticklabels(c_label, fontsize = font_size)  
        cbar.set_label(c_label_text, fontsize = font_size)   
   
def compute_streamlines(ax, v, points, vel_max, fig, font_size, i):
    print("Computing streamlines...")

    # --- Forward ---
    for i in range(len(points)):
        found = True
        j = 0

        px = points[i][0]
        py = points[i][1]

        trajectory_x = []
        trajectory_y = []
        vel_mag      = []

        while (found == True and j < 1000):
            v1 = v(Point(px,py))
            dt = 1e2/sqrt(v1[0]**2 + v1[1]**2)
            
            try:    
                v2 = v(Point(px + v1[0]*dt/2.0, py + v1[1]*dt/2.0))
                v3 = v(Point(px + v2[0]*dt/2.0, py + v2[1]*dt/2.0))
                v4 = v(Point(px + v3[0]*dt,     py + v3[1]*dt))
                
                px += (v1[0] + 2.0*v2[0] + 2.0*v3[0] + v4[0])/6.0*dt
                py += (v1[1] + 2.0*v2[1] + 2.0*v3[1] + v4[1])/6.0*dt

                trajectory_x.append(px)
                trajectory_y.append(py)
                vel_mag.append(sqrt(v(Point(px, py))[0]**2 + v(Point(px, py))[1]**2))
            
            except:
                found = False

            j += 1
        # ax.plot(trajectory_x, trajectory_y, c = "white", linewidth = "0.5", alpha=0.4)

        tx = np.array(trajectory_x)
        ty = np.array(trajectory_y)
        vv = np.array(vel_mag)

        pp = np.array([tx, ty]).T.reshape(-1, 1, 2)
        segments = np.concatenate([pp[:-1], pp[1:]], axis=1)

        norm = plt.Normalize(0, vel_max)
        lc1 = LineCollection(segments, cmap="GnBu", norm=norm)
        lc1.set_array(vv)
        lc1.set_linewidth(0.5)
        ax.add_collection(lc1)

    # --- Backward ---
    dtm = - dt
    for i in range(len(points)):
        found = True
        j = 0

        px = points[i][0]
        py = points[i][1]

        trajectory_x = []
        trajectory_y = []
        vel_mag      = []

        while (found == True and j < 1000):
            v1 = v(Point(px,py))
            dtm = -1e2/sqrt(v1[0]**2 + v1[1]**2)

            try:    
                v2 = v(Point(px + v1[0]*dtm/2.0, py + v1[1]*dtm/2.0))
                v3 = v(Point(px + v2[0]*dtm/2.0, py + v2[1]*dtm/2.0))
                v4 = v(Point(px + v3[0]*dtm,     py + v3[1]*dtm))
                
                px += (v1[0] + 2.0*v2[0] + 2.0*v3[0] + v4[0])/6.0*dtm
                py += (v1[1] + 2.0*v2[1] + 2.0*v3[1] + v4[1])/6.0*dtm

                trajectory_x.append(px)
                trajectory_y.append(py)
                vel_mag.append(sqrt(v(Point(px, py))[0]**2 + v(Point(px, py))[1]**2))
            
            except:
                found = False

            j += 1
        # ax.plot(trajectory_x, trajectory_y, c = "white", linewidth = "0.5", alpha=0.4)

        tx = np.array(trajectory_x)
        ty = np.array(trajectory_y)
        vv = np.array(vel_mag)

        pp = np.array([tx, ty]).T.reshape(-1, 1, 2)
        segments = np.concatenate([pp[:-1], pp[1:]], axis=1)

        norm = plt.Normalize(0, vel_max)
        lc1 = LineCollection(segments, cmap="GnBu", norm=norm)
        lc1.set_array(vv)
        lc1.set_linewidth(0.5)
        ax.add_collection(lc1)

    fig.subplots_adjust(right=0.9)
    cax = ax.inset_axes([1.3, 0, 0.03, 1])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="GnBu"), cax=cax, orientation='vertical')

    cbar.set_ticks(ticks=[0,  0.5*vel_max,  1*vel_max], labels=['0', "", "max"], fontsize = font_size)
    cbar.set_label("Velocity magnitude", fontsize = font_size) 

def scalebar(axis, length, sbar_length, corner_offset, fs):

    # Scale bar extent
    xx = [length - sbar_length - corner_offset, length - corner_offset]
    
    # Scale bar
    sbar = mpatches.Rectangle((xx[0], corner_offset + sbar_length/2*0.8), sbar_length, sbar_length/15.0, facecolor='black', zorder = 11)

    # Surrounding rectangle
    rect = mpatches.Rectangle((xx[0]*0.99, corner_offset), sbar_length + xx[0]*0.02, sbar_length/2*1.1, facecolor='white', zorder = 10)

    axis.add_patch(sbar)
    axis.add_patch(rect)
    axis.text((xx[0] + xx[1])/2.0, corner_offset +  sbar_length/6, str(int(sbar_length/1e3)) + " km", fontsize=fs, zorder = 11, horizontalalignment='center')

def plot_trajectory(ax, name, ID, sim_step, tail, tail_length):
    infile = open("data_" + name + "/tracers/trajectories/tracer_" + str(ID[0]) + ".dat", "r") 
    lines = infile.readlines() 

    xx = []
    yy = []
    alphas = []

    fade = colors.to_rgb(ID[1]) + (0.0,)
    mycolor = colors.LinearSegmentedColormap.from_list('my',[fade, ID[1]])

    for line in lines:
        sline = line.split("\t")
        step = float(sline[1])

        xx.append(float(sline[2]))
        yy.append(float(sline[3])) 
        if (tail == "comet"):
            alphas.append(max(0.0, (step - (sim_step - tail_length))/tail_length))

        if (step == sim_step):
            break
    
    if (tail == "comet"):
        points = np.vstack((xx, yy)).T.reshape(-1, 1, 2)
        segments = np.hstack((points[:-1], points[1:]))

        lc = LineCollection(segments, array = alphas, cmap = mycolor, lw = 3)
        line = ax.add_collection(lc)

    elif (tail == "full"):
        ax.plot(xx, yy, color = ID[1], linewidth = 2, zorder = 9) 
    
    else:
        pass

    ax.scatter(xx[-1], yy[-1],  color = ID[1], marker="o", edgecolor = "black", linewidths = 1, zorder = 10) 


def load_comp_tracers(ax, name, i, data):
    print("Plotting composition tracers...")
    xx = []
    yy = []

    for j in range(len(data)):
        xx.append([])
        yy.append([])

    infile = open("data_" + name + "/tracers/step_" + str(i) + ".dat", "r") 
    lines = infile.readlines() 

    header = True
    for line in lines:
        sline = line.split("\t")

        if (header == True):
            comp_idx = sline.index("composition")
            header = False
            continue

        for j in range(len(data)):
            if (float(sline[comp_idx + data[j][0]]) == 1.0):
                xx[j].append(float(sline[0]))
                yy[j].append(float(sline[1])) 

    for j in range(len(data)):
        ax.scatter(xx[j], yy[j],  c =  data[j][1], marker=".", s = 0.1, alpha=data[j][2]) 

def load_melt_tracers(ax, name, labels, fig, font_size, i):
    print("Plotting melt tracers...")
    xx = []
    yy = []
    mm = []

    m_label, m_label_text, m_label_val = labels[3][0], labels[3][1], labels[3][2]

    infile = open("data_"+name+"/tracers/step_"+str(i)+".dat", "r") 
    lines = infile.readlines() 

    header = True
    for line in lines:
        sline = line.split("\t")

        if (header==True):
            melt_idx = sline.index("xm (-)")
            header = False
            continue
        
        if (float(sline[melt_idx]) > 0.0):
            xx.append(float(sline[0]))
            yy.append(float(sline[1]))
            mm.append(float(sline[2])) 

    melt = ax.scatter(xx, yy,  c = mm, marker=",", s = 0.25, vmin=m_label_val[0], vmax=m_label_val[2], cmap="Blues") 

    # Colorbar 
    fig.subplots_adjust(right=0.9)
    cax = ax.inset_axes([1.3, 0, 0.03, 1])
    cbar = fig.colorbar(melt, cax=cax, orientation='vertical', ticks = m_label_val)

    cbar.ax.set_yticklabels(m_label, fontsize = font_size)  
    cbar.set_label(m_label_text, fontsize = font_size) 

def extract_cell_data(mesh, function, index):
    for j in range(mesh.num_cells()):                
            if (Cell(mesh,j).index() == index):
                centroid = Cell(mesh, j).midpoint()
                value = function(Point(centroid.x(), centroid.y()))   
    return value

def plot_time(ax, directories, column, labelname, colors, styles):
    for j in range(len(directories)):

        try:
            infile = open(directories[j], "r") 
            lines = infile.readlines() 

            time = []
            data_y = []

            header = True
            for line in lines:
                sline = line.split("\t\t")

                if (header == True):
                    header = False
                    continue
                else:
                    time.append(float(sline[0]))
                    data_y.append(float(sline[column])*1e3)
                    
            ax.plot(time, data_y, c = colors[j], lw = 2, linestyle = styles[j], label = labelname[j])

        except:
            print("Directory " + directories[j] + " not found.")
