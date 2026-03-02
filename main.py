# --- Python modules ---
from dolfin import *
import numpy as np
import time

# --- MilleFEuiIle modules ---
from m_mesh import *
from m_check import *
from m_incompatibility import *
from m_elements import *
from m_equations import *
from m_filenames import *
from m_tracers import *

from m_melting import *
from m_postproc import *
from m_timestep import *
from m_rheology import *
from m_constants import *
from m_parameters_docs import *
from m_interpolation import *
from m_material_properties import *
from m_boundary_conditions import *

def run_code():

    comm = MPI.comm_world
    rank = MPI.rank(comm)
    size = MPI.size(comm)

    code_start = time.time()
    code_now_k = time.time()

    if (rank == 0):
        print(" ", flush=True)
        print("\t"+"~ M i l l e F E u i I l e ~", flush=True)
        print("\t"+"___________     ___________", flush=True)
        print("\t"+"(__  __)::::\\__/:::::::::::", flush=True)
        print("\t"+"   )) (__ __)\\/(    )(    )", flush=True)
        print("\t"+"  ((     ((   :::::::::::::", flush=True)
        print("\t"+"   ))     ))  (    )(     )", flush=True)
        print("\t"+":::::::::::::::::::::::::::", flush=True)
        print(" ", flush=True)

    # --- Check of the basic string inputs ---
    check_input_parameters()
    check_compatibility()

    MPI.barrier(comm)

    t = Constant(0.0)
    step = 0

    # --- Import classes ---
    MeshClass       = MeshModule()
    ElemClass       = Elements(MeshClass)
    FilesClass      = SaveFiles(MeshClass, ElemClass)
    TracersClass    = Tracers(MeshClass, ElemClass, FilesClass)
    MeltingClass    = Melting(MeshClass, ElemClass, TracersClass)

    EqClass         = Equations(MeshClass, ElemClass, TracersClass, MeltingClass, FilesClass, t)

    if (termination_condition[0] != "initial_condition"):
        if (termination_condition[0] == "time"):
            t_end = termination_condition[1]
        else:
            step_end = termination_condition[1]

        if (output_frequency[0] == "time"):
            time_output  = np.linspace(0, t_end, int(t_end/output_frequency[1] + 1))
        else:
            time_output  = None

    step_output = 0
        
    if (reload == True):
        if (solve_energy_problem == True):
            # --- Load temperature ---
            FilesClass.Load_HDF5(t, EqClass.dt, "temperature", EqClass.Temp_k)
            EqClass.Temp.assign(EqClass.Temp_k)
        
        if (BC_Stokes_problem[0][0] == "free_surface"):
            # --- Load topography and assign BCs for the mesh displacement ---
            FilesClass.Load_HDF5(t, EqClass.dt, "topography_top", EqClass.h_top)
            EqClass.h2_top.assign(EqClass.h_top)

        if (BC_Stokes_problem[1][0] == "free_surface"):
            FilesClass.Load_HDF5(t, EqClass.dt, "topography_bottom", EqClass.h_bot)
            EqClass.h2_bot.assign(EqClass.h_bot)

        if (BC_Stokes_problem[0][0] == "free_surface" or BC_Stokes_problem[1][0] == "free_surface"):  
            # --- Move the mesh ---
            EqClass.solver_surf_move.solve()
            EqClass.u_mesh.assign(project(EqClass.h2*EqClass.e_z, ElemClass.vCG1))
            EqClass.h_top_k.assign(EqClass.h_top)
            EqClass.h_bot_k.assign(EqClass.h_bot)

            MeshClass.move_mesh(ElemClass.u_mesh)

        if (TracersClass.use_tracers == True):
            TracersClass.load_tracers()

            # --- Interpolate tracer-carried functions ---
            for i in range(len(materials)):
                composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])
            
            tracer_count_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells,  EqClass.number_of_tracers)

            if (internal_melting == True):
                melt_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 0, EqClass.xm_k)
            
            if (plasticity == True):
                scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)
        
            # if (elasticity == True):
            #     EqClass.update_stress()

    else:
        EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))

        if (initial_topography == True):
            EqClass.h_top.assign(project(h_top_ini, ElemClass.sCG1))
            EqClass.h_bot.assign(project(h_bot_ini, ElemClass.sCG1))

            EqClass.h2_top.assign(EqClass.h_top)
            EqClass.h2_bot.assign(EqClass.h_bot)

            EqClass.solver_surf_move.solve()
            EqClass.u_mesh.assign(project(EqClass.h2*EqClass.e_z, ElemClass.vCG1))
            EqClass.h_top_k.assign(EqClass.h_top)
            EqClass.h_bot_k.assign(EqClass.h_bot)

            MeshClass.move_mesh(ElemClass.u_mesh)

        if (TracersClass.use_tracers == True):
            # TracersClass.introduce_tracers_unstructured()
            TracersClass.tracer_count_interpolation()
            scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)

        # Define the temperature even though it will not be solved
        EqClass.Temp.assign(project(Expression("Tb - x[1]/h*(Tb-Ts)", Tb=BC_heat_transfer[1][1], Ts=BC_heat_transfer[0][1], h=height, degree=1), ElemClass.sCG2))

        if (solve_energy_problem == True):             
            
            if (init_cond_profile == True):
                EqClass.solve_initial_heat_equation()
            
            if (cos_perturbation == True):
                EqClass.thermal_perturbation()

    # --- Interpolate composition ---
    for i in range(len(materials)):
        composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])

    # --- Save initial condition ---
    FilesClass.Save_Paraview_Ini()

    if (termination_condition[0] == "initial_condition"):
        if (rank == 0):
            print("Exiting after initial condition.")
        exit()
        
    # --- Main time loop ---
    while ((termination_condition[0] == "time" and float(t) < t_end)
        or (termination_condition[0] == "step" and step < step_end)):
        
        # --- Check cache ---
        if (monitor_cache == True and  rank == 0):
            print("\nCache control", flush=True)
            os.system("ls /home/kihoulou/.cache/dijitso/lib/ | wc -l\n", flush=True)

        # --- Solve Stokes problem ---
        EqClass.solve_Stokes_problem(step, FilesClass)

        # --- Update time and step ---
        t.assign(float(t + EqClass.dt))
        step += 1

        MPI.barrier(comm)
        if (rank == 0):
            print("\n----------------------------------------------", flush=True)
            print("\tStep:     ", '{:d}'.format(step), flush=True)
            print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string, flush=True)
            print("\tTime step:", '{:.3e}'.format(float(EqClass.dt/time_units)), time_units_string, flush=True)
            print("----------------------------------------------\n", flush=True)

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

                # if (step == 1):
                #     TracersClass.introduce_tracers_unstructured()
                # TracersClass.add_tracers_unstructured()

                TracersClass.tracer_count_interpolation()
                

            # --- Interpolate tracer-carried functions ---
            for i in range(len(materials)):
                composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i])
            
            tracer_count_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells,  EqClass.number_of_tracers)

            scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)
        
        # --- Update stress ---
        if (elasticity == True):
            EqClass.update_stress()
            EqClass.rotate_and_interpolate_stress()

        # --- Postprocessing ---
        EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
        EqClass.q_top.assign(assemble(dot(-k(EqClass.Temp, EqClass.composition)*nabla_grad(EqClass.Temp), EqClass.normal)*MeshClass.ds(1)) / EqClass.top_length)
        EqClass.log10_visc.assign(project(ln(EqClass.visc)/ln(10.0), ElemClass.sDG0))
        EqClass.cohesion.assign(project(cohesion(EqClass.plastic_strain), ElemClass.sDG0))

        # --- Check whether to save results? ---> If yes, save them.
        step_output, output_now = Output_Timing(step, step_output, t, time_output)
        if (output_now == True):
            if (TracersClass.use_tracers == True):
                # TracersClass.tracer_count_interpolation()
                TracersClass.save_tracers(step_output)

            TracersClass.rank_interpolation()
            FilesClass.Save_Paraview(t)
            FilesClass.Save_HDF5(step_output, step, EqClass.dt, t)

        code_now        = time.time()
        total_time      = (code_now - code_start)/3600.0
        timestep_time   = (code_now - code_now_k)
        code_now_k      = time.time()

        FilesClass.write_statistic(t, step, stat_output,\
                                q_cond_top  = EqClass.q_cond_top,\
                                q_top       = EqClass.q_top,\
                                v           = EqClass.v_k,\
                                avg_h_bot   = EqClass.h_bot_aver,\
                                h_top_max   = EqClass.h_top,\
                                time        = total_time,\
                                timestep    = timestep_time)
        
        # --- Frozen ice shell termination criterium ---
        # if (float(EqClass.thickness_now) > 167e3):
        #     if (rank == 0):
        #         print("\n----------------------------------------------")
        #         print("\tSubsurface ocean completely frozen.")
        #         print("\tStep:     ", '{:d}'.format(step))
        #         print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string)
        #         print("----------------------------------------------\n")
        #     break

        # --- Equilibrated surface heat flux termination criterium ---
        # if (float(t) > 3*Myr and step % 10 == 0):
        #     infile = open("data_" + name + "/statistics.dat", "r") 
        #     lines = infile.readlines() 

        #     read_time = []
        #     read_qtop = []

        #     header = True
        #     for line in lines:
        #         sline = line.split("\t\t")

        #         if (header==True):
        #             header = False
        #             continue

        #         read_time.append(float(sline[0]))
        #         read_qtop.append(float(sline[2])) 

        #     qtop1, qtop2, qtop3 = [], [], []
        #     for i in range(0, len(read_time)-1):
        #         if (read_time[i] > read_time[-1] - 1):
        #             qtop1.append(read_qtop[i])

        #         if (read_time[i] > read_time[-1] - 2):
        #             qtop2.append(read_qtop[i])

        #         if (read_time[i] > read_time[-1] - 3):
        #             qtop3.append(read_qtop[i])

        #     qtop1_aver = sum(qtop1) / len(qtop1)
        #     qtop2_aver = sum(qtop2) / len(qtop2)
        #     qtop3_aver = sum(qtop3) / len(qtop3)
            
        #     tol_q = 1e-3
        #     tol_v = 1e-2
        #     if ((1.0 + tol_q)*qtop3_aver > qtop1_aver > (1.0 - tol_q)*qtop3_aver and (1.0 + tol_q)*qtop3_aver > qtop2_aver > (1.0 - tol_q)*qtop3_aver):
        #         equilibrated = True
        #     else:
        #         equilibrated = False

        #     if (equilibrated == True):
        #         if (rank == 0):
        #             print("\n----------------------------------------------")
        #             print("\tSteady state reached.")
        #             print("\tStep:     ", '{:d}'.format(step))
        #             print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string)
        #             print("----------------------------------------------\n")
        #         break
        
    FilesClass.Save_Paraview(t)
    FilesClass.Save_HDF5(step_output, step, EqClass.dt, t)

    if (rank == 0):
        print("\n----------------------------------------------")
        print("\tTermination condition satisfied.")
        print("\tStep:     ", '{:d}'.format(step))
        print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string)
        print("----------------------------------------------\n")

if __name__ == "__main__":
    run_code()