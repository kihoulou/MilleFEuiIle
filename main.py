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
from m_parameters import *
from m_interpolation import *
from m_material_properties import *
from m_boudnary_conditions import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

code_start = time.time()
code_now_k = time.time()

if (rank == 0):
    print(" ")
    print("\t"+"~ M i l l e F E u i I l e ~")
    print("\t"+"___________     ___________")
    print("\t"+"(__  __)::::\\__/:::::::::::")
    print("\t"+"   )) (__ __)\\/(    )(    )")
    print("\t"+"  ((     ((   :::::::::::::")
    print("\t"+"   ))     ))  (    )(     )")
    print("\t"+":::::::::::::::::::::::::::")
    print(" ")

# --- Check of the basic string inputs ---
if (rank == 0):
    check_input_parameters()

MPI.barrier(comm)

# --- Import classes ---
MeshClass       = MeshModule()
ElemClass       = Elements(MeshClass)
TracersClass    = Tracers(MeshClass, ElemClass)
MeltingClass    = Melting(MeshClass, ElemClass, TracersClass)

FilesClass      = SaveFiles(MeshClass, ElemClass)
EqClass         = Equations(MeshClass, ElemClass, TracersClass, MeltingClass)

# --- Save the source code ---
if (rank == 0):
    os.system("cp  main.py data_" + name + "/source_code")
    os.system("cp  m_*.py data_" + name + "/source_code")

t = Constant(0.0)
step = 0

if (output_frequency[1] == "time"):
    time_output  = np.linspace(0, t_end, int(t_end/output_frequency[1] + 1))
else:
    time_output  = None
step_output = 0

        
if (reload_HDF5 == True):
    FilesClass.Load_HDF5(t, EqClass.dt, EqClass.Temp_k)
    EqClass.Temp.assign(EqClass.Temp_k)

else:
    EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
    if (solve_energy_problem == True):
        EqClass.Temp.assign(project(Expression("Tb - x[1]/h*(Tb-Ts)", Tb=T_bot, Ts=T_top, h=height, degree=1), ElemClass.sCG2))             
        if (init_cond_profile == True):
            EqClass.solve_initial_heat_equation()
        
        if (cos_perturbation == True):
            EqClass.thermal_perturbation()

# --- Interpolate composition ---
for i in range(len(materials)):
    composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])


# --- Guess the phase boundary velocity a priori ---
# if (BC_vel_bot == "free_surface" and phase_transition == True):

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
    EqClass.solve_Stokes_problem(step, FilesClass)

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
        EqClass.solve_topography_evolution()

    # --- Advect tracers and move mesh ---
    MPI.barrier(comm)
    # Advection part I
    if (TracersClass.use_tracers == True):
        TracersClass.advect_tracers(EqClass.v_k, EqClass.v_mesh, EqClass.dt)
        
    # Mesh displacement
    if (MeshClass.moving_mesh == True):
        MeshClass.move_mesh(ElemClass.u_mesh)

    # Advection part II
    if (TracersClass.use_tracers == True):
        TracersClass.find_tracers(EqClass.v_k, EqClass.dt)        
        TracersClass.delete_and_find()

        if (TracersClass.only_melt_tracers == False):
            TracersClass.add_tracers(t, step)

        # --- Interpolate tracer-carried functions ---
        for i in range(len(materials)):
            composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])
    
    # --- Update stress ---
    if (elasticity == True):
        EqClass.update_stress()

    # --- Postprocessing ---
    EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
    EqClass.q_top.assign(assemble(dot(-k(EqClass.Temp, EqClass.composition)*nabla_grad(EqClass.Temp), EqClass.normal)*MeshClass.ds(1)) / EqClass.top_length)
    EqClass.log10_visc.assign(project(ln(EqClass.visc)/ln(10.0), ElemClass.sDG0))

    # --- Check whether to save results? ---> If yes, save them.
    step_output, output_now = Output_Timing(step, step_output, t, time_output)
    if (output_now == True):
        TracersClass.rank_interpolation()
        FilesClass.Save_Paraview(t)
        FilesClass.Save_HDF5(step_output, step, EqClass.dt, t)

        if (TracersClass.use_tracers == True):
            TracersClass.tracer_count_interpolation()
            if (save_tracers == True):
                TracersClass.save_tracers(step)

    code_now        = time.time()
    total_time      = (code_now - code_start)/3600.0
    timestep_time   = (code_now - code_now_k)/3600.0
    code_now_k      = time.time()

    FilesClass.write_statistic(t, step, stat_output,\
                            q_cond_top  = EqClass.q_cond_top,\
                            q_top       = EqClass.q_top,\
                            v           = EqClass.v_k,\
                            avg_h_bot   = EqClass.h_bot_aver,\
                            time        = total_time,\
                            timestep    = timestep_time)
    
FilesClass.Save_Paraview(t)
FilesClass.Save_HDF5(step_output, step, EqClass.dt, t)