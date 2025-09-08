from dolfin import *
import matplotlib.tri as tri
import matplotlib.patches as mpatches

comm = MPI.comm_world
rank = MPI.rank(comm)
size= MPI.size(comm)

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def initialize(name, i):
    mesh = Mesh()

    mesh_file = HDF5File(comm, "data_"+name+"/HDF5/meshes/mesh_" + str(i) +".h5", "r")
    mesh_file.read(mesh, "/mesh", True)
    mesh_file.close()

    sDG0_elem = FiniteElement("DG", mesh.ufl_cell(), 0)
    sDG0 = FunctionSpace(mesh, sDG0_elem)

    sCG2_elem = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
    sCG2 = FunctionSpace(mesh, sCG2_elem)

    vCG1_elem = VectorElement("Lagrange", mesh.ufl_cell(), 1)
    vCG1 = FunctionSpace(mesh, vCG1_elem)

    Temp = Function(sCG2)
    v = Function(vCG1)
    visc = Function(sDG0)
    visc_log = Function(sDG0)

    return Temp, v, visc, mesh, visc_log, sDG0

def load(i, Temp, name):
    Temp_dataset = "/temperature/vector_%d"% (i)
    v_dataset = "/velocity/vector_%d"% (i)
    visc_dataset = "/viscosity/vector_%d"% (i)
    
    files_hdf5_in={}
    files_hdf5_in['data'] = HDF5File(comm, "data_"+name+"/HDF5/data.h5", "r")

    files_hdf5_in['data'].read(Temp, Temp_dataset)
    # files_hdf5_in['data'].read(v, v_dataset)
    # files_hdf5_in['data'].read(visc, visc_dataset)

    return Temp #, v, visc

def plot_scalar(ax, fig, function, mesh, labels, x_range, z_range, font_size):
    
    ax.set_aspect("equal")
    mesh_function = function.function_space().mesh()
    C = function.compute_vertex_values(mesh)

    values = []
    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()
        values.append(function(Point(centroid.x(), centroid.y())))
    
    # tpc = ax.tripcolor(mesh2triang(mesh_function), values, shading='flat', vmin = 14, vmax = 23, cmap = "jet_r")
    
    tpc = ax.tripcolor(mesh2triang(mesh_function), C, shading='gouraud',  vmin = 100, vmax = 273, cmap = "coolwarm")

    # Axes range
    ax.set_ylim(z_range[0],z_range[1])
    ax.set_xlim(x_range[0],x_range[1])

    x_label, x_label_text, x_label_val = labels[0][0], labels[0][1], labels[0][2]
    z_label, z_label_text, z_label_val = labels[1][0], labels[1][1], labels[1][2]
    c_label, c_label_text, c_label_val = labels[2][0], labels[2][1], labels[2][2]

    # x-label settings
    ax.set_xlabel(x_label_text, fontsize = font_size)
    ax.set_xticks(x_label_val)
    ax.set_xticklabels(x_label, fontsize = font_size)

    # y-label settings
    ax.set_ylabel(z_label_text, fontsize = font_size)
    ax.set_yticks(z_label_val)
    ax.set_yticklabels(z_label, fontsize = font_size)

    ax.tick_params(pad=10) 

    # Colorbar 
    fig.subplots_adjust(right=0.9)
    cax = ax.inset_axes([1.05, 0, 0.03, 1])
    cbar = fig.colorbar(tpc, cax=cax, orientation='vertical', ticks = c_label_val)

    cbar.ax.set_yticklabels(c_label, fontsize = font_size)  
    cbar.set_label(c_label_text, fontsize = font_size)   
   
def compute_streamlines(ax, v, points):
    print("Computing streamlines...")

    # --- Forward ---
    for i in range(len(points)):
        found = True
        j = 0

        px = points[i][0]
        py = points[i][1]

        trajectory_x = []
        trajectory_y = []

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
            
            except:
                found = False

            j += 1
        ax.plot(trajectory_x, trajectory_y, c = "white", linewidth = "0.5", alpha=0.4)

     # --- Backward ---
    dtm = - dt
    for i in range(len(points)):
        found = True
        j = 0

        px = points[i][0]
        py = points[i][1]

        trajectory_x = []
        trajectory_y = []

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
            
            except:
                found = False

            j += 1
        ax.plot(trajectory_x, trajectory_y, c = "white", linewidth = "0.5", alpha=0.4)

def scalebar(axis, length, y_coord, fs):
    right_off = 2e3
    sbar_length = 10e3
    xx = [length - sbar_length - right_off, length - right_off]
    
    # Scale bar
    sbar = mpatches.Rectangle((xx[0], y_coord + 4.5e3), sbar_length, 5e2, facecolor='black', zorder = 11)

    # Surrounding rectangle
    rect = mpatches.Rectangle((xx[0] - 1e3, y_coord + 1e3), sbar_length + 2e3, 5e3, facecolor='white', zorder = 10)

    axis.add_patch(sbar)
    axis.add_patch(rect)
    axis.text((xx[0] + xx[1])/2.0, y_coord + 2e3, "10 km", fontsize=fs, zorder = 11, horizontalalignment='center')

def load_tracers(ax, name, i):
    print("Plotting tracers...")
    xx_0 = []
    yy_0 = []

    xx_1 = []
    yy_1 = []

    infile = open("data_"+name+"/tracers/step_"+str(i)+".dat", "r") 
    lines = infile.readlines() 

    header = True
    for line in lines:
        sline = line.split("\t")

        if (header==True):
            header = False
            continue
        
        if (float(sline[2]) == 1.0 and float(sline[1]) > 24e3):
            xx_0.append(float(sline[0]))
            yy_0.append(float(sline[1])) 
        
        if (float(sline[3]) == 1.0):
            xx_1.append(float(sline[0]))
            yy_1.append(float(sline[1])) 

    ax.scatter(xx_0, yy_0,  c = "lightgrey", marker=",", s = 1.0, alpha=1) 
    ax.scatter(xx_1, yy_1,  c = "black", marker=",", s = 1.0, alpha=1) 

def load_melt_tracers(ax, name, labels, fig, font_size, i):
    m_label, m_label_text, m_label_val = labels[3][0], labels[3][1], labels[3][2]

    print("Plotting tracers...")
    xx = []
    yy = []
    mm = []

    infile = open("data_"+name+"/tracers/step_"+str(i)+".dat", "r") 
    lines = infile.readlines() 

    header = True
    for line in lines:
        sline = line.split("\t")

        if (header==True):
            header = False
            continue
        
        if (float(sline[2]) > 0.0):
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
