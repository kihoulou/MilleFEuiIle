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
        or (output_frequency[0] == "steps" and step % output_frequency[1] == 0)\
        or (output_frequency[0] == "time" and float(t) >= time_output[step_output])):
        value = True
        step_output += 1
        
    else:
        value = False

    return step_output, value

def time_step(mesh, v, v_mesh, H_max, composition, Temp, unit_scalar, t):
    

    if (time_step_strategy == "constant"):
        return dt_const
    
    if (time_step_strategy == "domain"):
        # --- Determine minimum elememt size in the mesh ---
        x_min_ranks = length
        for f in facets(mesh):
                for e in edges(f):
                    if (e.length()< x_min_ranks): x_min_ranks = e.length()
        
        x_min = MPI.min(comm, x_min_ranks) 

        # --- Compute the maximum speed in the domain ---
        v_max = MPI.max(mesh.mpi_comm(), np.abs(v.vector().get_local()).max()) 
        
        dt_list = []

        # --- Convective time step ---
        dt_conv = cfl*x_min/v_max
        # dt_list.append(dt_conv)

        # --- Time-scaled time step ---
        if (scaled_time_step == True):
            # Useful for problems with logaritmic time evolution, such as benchmark Patocka et al. (2017)
            dt_time = float(t)/time_step_scaling
            dt_list.append(dt_time)

        # --- Fixed maximum time step ---
        # if (maximum_time_step == True):
        #     dt_list.append(dt_max)

        # --- Conductive time step ---
        if (solve_energy_problem == True):
            temp_aver = assemble(Temp*dx)/assemble(unit_scalar*dx)
            dt_cond = cfl*x_min**2*rho_s*cp(temp_aver, composition)/k(temp_aver, composition)
            dt_list.append(dt_cond)

        # --- Mesh displacement time step ---
        if (BC_Stokes_problem[0][0] == "free_surface" or BC_Stokes_problem[1][0] == "free_surface"):
            v_max_mesh = MPI.max(mesh.mpi_comm(), np.abs(v_mesh.vector().get_local()).max()) 
            dt_mesh = cfl*(height/1e3)/(v_max_mesh + 1e-15)
            dt_list.append(dt_mesh)
        
        # --- Internal heating time step ---
        if (tidal_dissipation == True):
            dt_H = cfl*rho_s*cp(temp_aver, composition)*dT_max/H_max
            dt_list.append(dt_H)

        return min(dt_list)