from dolfin import *
import numpy as np

from m_constants import *
from m_parameters import *
from m_material_properties import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

def Output_Timing(step, step_output, t, time_output):
    """Determines whether to save the data in the current time step.

    :param step: order number of the current step in the time loop (increases by 1 with every time step)
    :param step_output: order number of the output (increases by 1 with every output)
    :param t: time
    :param time_output: list of time values at which the data will be saved.

    :returns: ``True`` or ``False`` and an updated value of the ``step_output`` parameter.

    """
    
    if (step == 1\
        or (output_frequency[0] == "steps" and step % output_frequency[1] == 0)\
        or (output_frequency[0] == "time" and float(t) >= time_output[step_output])):
        value = True
        step_output += 1
        
    else:
        value = False

    return step_output, value

def time_step(mesh, v, v_mesh, H_max, composition):
    """Determines the length of the next time step :math:`\\Delta t`\ .

    :param mesh: computational mesh
    :param v: material velocity (:math:`v`\ )
    :param v_mesh: mesh displacement velocity (:math:`v_{mesh}`\ )
    :param H_max: peak tidal heating
    :param composition: material composition

    :returns: ahoj
    """
    

    if (time_step_strategy == "constant"):
        """Returns a constant time step. """
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
        dt_list.append(dt_conv)

        # --- Conductive time step ---
        if (solve_energy_problem == True):
            dt_cond = cfl*x_min**2*rho_s*cp((T_bot+T_top)/2.0, composition)/k((T_bot+T_top)/2.0 , composition)
            dt_list.append(dt_cond)

        # --- Mesh displacement time step ---
        if (BC_vel_bot == "free_surface" or BC_vel_top == "free_surface"):
            v_max_mesh = MPI.max(mesh.mpi_comm(), np.abs(v_mesh.vector().get_local()).max()) 
            dt_mesh = cfl*(height/1e3)/(v_max_mesh + 1e-15)
            dt_list.append(dt_mesh)
        
        # --- Internal heating time step ---
        if (tidal_dissipation == True):
            dt_H = cfl*rho_s*cp(T_bot, composition)*dT_max/H_max
            dt_list.append(dt_H)

        return min(dt_list)