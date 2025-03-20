# --- Python modules ---
from dolfin import *
import numpy as np
import time

# --- MilleFEuiIle modules ---
from m_mesh import *
from m_check import *
from m_elements import *
from m_equations import *
from m_filenames import *
from m_tracers import *

from m_melting import *
from m_postproc import *
from m_timestep import *
from m_rheology import *
from m_constants import *

# from m_advection import *
from m_parameters import *
from m_interpolation import *
from m_material_properties import *
from m_boudnary_conditions import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

code_start = time.time()
code_now_k = time.time() #Ahoj

if (rank == 0):
    print(" ")
    print("\t"+"~ M i l l e F E u i I l e ~")
    print("\t"+"___________     ___________")
    print("\t"+"(__  __)::::\\__/:::::::::::")
    print("\t"+"   )) (__ __)\\/(    )(    )")
    print("\t"+"  ((     ((   :::::::::::::")
    print("\t"+"   ))     ))  (    )(     )")
    print("\t"+":::::::::::::::::::::::::::xx")
    print(" ")

# --- Check of the basic string inputs ---
if (rank == 0):
    Check_Input_Parameters()

# --- Import classes ---
MeshClass       = MeshModule()
ElemClass       = Elements(MeshClass)
EqClass         = Equations(MeshClass, ElemClass)
FilesClass      = SaveFiles(MeshClass, EqClass)

t = Constant(0.0)
step = 0

time_output  = np.linspace(0, t_end, int(t_end/every_t + 1))
step_output = 0
        
if (reloading_HDF5 == True):
    FilesClass.Load_HDF5(t, EqClass.dt, EqClass.Temp_k)
    EqClass.Temp.assign(EqClass.Temp_k)

else:
    EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
    if (solve_energy_problem == True):
        EqClass.solver_energy_ini.solve()
        EqClass.q_cond_top.assign(assemble(dot(-k(EqClass.Temp, EqClass.composition)*nabla_grad(EqClass.Temp), EqClass.normal)*MeshClass.ds(1)) / EqClass.top_length)

        if (cos_perturbation == True):
            EqClass.PerturbProfile()

# --- Import tracer class once the mesh is set and we know the temperature
TracersClass    = Tracers(MeshClass, EqClass)

# --- Interpolate composition ---
for i in range(len(materials)):
    composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])

# --- Main time loop structure ---

# 1/ solve Stokes problem
# 2/ update time step
# 3/ solve energy problem
# 4/ solve topography evolution
# 5/ advect tracers and move mesh
# 6/ update tracer-carried functions
# 7/ save data

while (float(t) < t_end):
    
    # --- Solve Stokes problem ---
    EqClass.solve_Stokes_problem()

    # --- Update time and step ---
    t.assign(float(t + EqClass.dt))
    step += 1

    MPI.barrier(comm)
    if (rank == 0):
        print(" ")
        print("----------------------------------------------")
        print("\tStep:     ", '{:d}'.format(step))
        print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string)
        print("\tTime step:", '{:.3e}'.format(float(EqClass.dt/time_units)), time_units_string)
        print("----------------------------------------------")
        print(" ")

    # --- Solve heat transfer equation ---
    if (solve_energy_problem == True):
        EqClass.solve_heat_equation()

    # --- Solve motion of the boundaries ---
    if (MeshClass.moving_mesh == True):
        EqClass.Solver_Topography_evolution()

    # --- Advect tracers and move mesh ---
    MPI.barrier(comm)
    if (TracersClass.use_tracers == True):
        TracersClass.advect_tracers(EqClass.v_k, EqClass.v_mesh, EqClass.dt)
        
        if (MeshClass.moving_mesh == True):
            MeshClass.move_mesh(EqClass.u_mesh)

        TracersClass.find_tracers(EqClass.v_k, EqClass.dt)        
        TracersClass.delete_and_find()
        TracersClass.add_tracers(t, step)

    # --- Interpolate tracer-carried functions ---
    for i in range(len(materials)):
        composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])
    
    # --- Check whether to save results? ---> If yes, save them.
    step_output, output_now = Output_Timing(step, step_output, t, time_output)
    if (output_now == True):
        TracersClass.tracer_count_interpolation()
        TracersClass.rank_interpolation()

        FilesClass.Save_Paraview(t)
        FilesClass.Save_HDF5(step_output, step, EqClass.dt, t)

        if (save_tracers == True):
            TracersClass.save_tracers(step)

        # --- Postprocessing 
    EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
    EqClass.q_top.assign(assemble(dot(-k(EqClass.Temp, EqClass.composition)*nabla_grad(EqClass.Temp), EqClass.normal)*MeshClass.ds(1)) / EqClass.top_length)

    code_now        = time.time()
    total_time      = (code_now - code_start)/3600.0
    timestep_time   = (code_now - code_now_k)/3600.0
    code_now_k      = time.time()

    
    FilesClass.write_statistic(t, step, stat_output,\
                            q_cond_top  = EqClass.q_cond_top,\
                            q_top       = EqClass.q_top,\
                            v           = EqClass.v_k,\
                            time        = total_time,\
                            timestep    = timestep_time)
    
FilesClass.Save_Paraview(t)
FilesClass.Save_HDF5(step_output)