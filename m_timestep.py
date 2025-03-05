from dolfin import *
import numpy as np

from m_constants import *
from m_parameters import *
from m_material_properties import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

def Output_Timing(step, step_output, t, time_output):
    if (step == 1\
        or (output_type == "steps" and step % every_n == 0)\
        or (output_type == "time" and float(t) >= time_output[step_output])):
        value = True
        step_output += 1
        
    else:
        value = False

    return step_output, value

def time_step_domain(mesh, v, H_max, composition):

    x_min_ranks = length
    for f in facets(mesh):
            for e in edges(f):
                if (e.length()< x_min_ranks): x_min_ranks = e.length()
    
    x_min = MPI.min(comm, x_min_ranks) 
    v_max = MPI.max(mesh.mpi_comm(), np.abs(v.vector().get_local()).max()) 

    dt_conv = cfl*x_min/v_max
    dt_cond = cfl*x_min**2*rho_s*cp((T_bot+T_top)/2.0, composition)/k((T_bot+T_top)/2.0 , composition)
    dt_H_max = cfl*rho_s*cp(T_bot, composition)*dT_max/H_max
    

    return dt_conv
    # if (solve_energy_problem == True):
    #     if (tidal_dissipation == True):
    #         return min(min(dt_conv, dt_cond), dt_H_max)
    #     else:
    #         return min(dt_conv, dt_cond)
    # else:
    #     return dt_conv
    
def time_step_cell(v, Temp, timestep):
    ranks = []
    dt_min = 10.0*Myr

    for j in range(mesh.num_cells()):
        cell_j = Cell(mesh,j).get_vertex_coordinates()
        centroid = Cell(mesh, j).midpoint()

        temp_cell = Temp(Point(centroid.x(),centroid.y()))
        v_cell = sqrt(v(Point(centroid.x(),centroid.y()))[0]**2 + v(Point(centroid.x(),centroid.y()))[1]**2)

        side_1 = sqrt((cell_j[0] - cell_j[2])**2 + (cell_j[1] - cell_j[3])**2)
        side_2 = sqrt((cell_j[2] - cell_j[4])**2 + (cell_j[3] - cell_j[5])**2)
        side_3 = sqrt((cell_j[0] - cell_j[4])**2 + (cell_j[1] - cell_j[5])**2)

        x_min = min(side_1,min(side_2, side_3))

        # --- Timestep for advection problems ---
        dt_conv = cfl*x_min/v_cell

        # --- Timestep for heat equation ---
        dt_cond = cfl*x_min**2*rho_s*cp(temp_cell)/k(temp_cell)

        min_timestep = min(dt_cond, dt_conv)
        ranks.append(min_timestep/kyr)

        if (min_timestep < dt_min):
            dt_min = min_timestep

    MPI.barrier(comm)
    dt_min = MPI.min(comm, dt_min)

    ranks = numpy.array(ranks)                                                            
    timestep.vector().set_local(ranks)
    
    return dt_min