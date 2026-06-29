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
from m_parameters import *
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
        print(" ", flush = True)
        print("\t"+"~ M i l l e F E u i I l e ~", flush = True)
        print("\t"+" ________________________", flush = True)
        print("\t"+"|__xxxxxxx_______________|", flush = True)
        print("\t"+"|  _______    ---------  |", flush = True)
        print("\t"+"| (__   __)  |	       | |	", flush = True)
        print("\t"+"|    \ /      ---------  |", flush = True)
        print("\t"+"|     \\\\     ::::::::::: |", flush = True)
        print("\t"+"|      \\\\     ---------  |", flush = True)
        print("\t"+"|      //    |         | |", flush = True)
        print("\t"+"|     / \     ---------  |", flush = True)
        print("\t"+"|__.-'   '--..___________|", flush = True)
        print("\t"+"|__.-'''''--..___________|", flush = True)
        print(" ", flush = True)

    # --- Check of the basic string inputs ---
    check_input_parameters()
    check_compatibility()

    MPI.barrier(comm)

    if (time_units_string == "-" or time_units_string == "s"):
        time_units = 1.0
    elif (time_units_string == "hr"): 
        time_units = hr
    elif (time_units_string == "day"): 
        time_units = day
    elif (time_units_string == "yr"): 
        time_units = yr
    elif (time_units_string == "kyr"): 
        time_units = kyr
    elif (time_units_string == "Myr"): 
        time_units = Myr   
    elif (time_units_string == "Gyr"): 
        time_units = Gyr

    # --- Determine whether to use tracers or not ---
    if (plasticity == True or elasticity == True or internal_melting == True or len(materials) > 0 or len(save_tracer_trajectory) > 0):
        use_tracers = True
    else:
        use_tracers = False

    # --- Determine whether the rheology is linear or not ---
    if (viscosity_type == "GK_2001" or plasticity == True):
        nonlinear_rheology = True
    else:
        nonlinear_rheology = False
        
    t = Constant(0.0)
    step = 0

    # --- Import classes ---
    MeshClass       = MeshModule()
    ElemClass       = Elements(MeshClass)
    FilesClass      = SaveFiles(MeshClass, ElemClass, use_tracers)
    TracersClass    = Tracers(MeshClass, ElemClass, FilesClass, use_tracers)
    MeltingClass    = Melting(MeshClass, ElemClass, TracersClass)

    EqClass         = Equations(MeshClass, ElemClass, TracersClass, MeltingClass, FilesClass, t, use_tracers, size)

    step_output = 0
    
    if (reload == True):
        if (solve_energy_problem == True):
            # --- Load temperature ---
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "temperature", EqClass.Temp)
            EqClass.Temp_k.assign(EqClass.Temp)
        
        if (BC_Stokes_problem[0][0] == "free_surface"):
            # --- Load topography and assign BCs for the mesh displacement ---
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "topography_top", EqClass.h_top)
            EqClass.h_top_k.assign(EqClass.h_top)

        if (BC_Stokes_problem[1][0] == "pressure"):
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "topography_bottom", EqClass.h_bot)
            EqClass.h_bot_k.assign(EqClass.h_bot)
        
        if (nonlinear_rheology == True):
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "velocity", EqClass.v_k)
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "viscosity", EqClass.visc)

        if (plasticity == True):
            FilesClass.load_HDF5(t, EqClass.dt, time_units, "pressure", EqClass.p_k)

        if (use_tracers == True):
            TracersClass.load_tracers()

            # --- Interpolate tracer-carried functions ---
            # - Composition -
            n_tracers_current = []
            n_tracers_orig = []
            for i in range(len(materials)):
                n_current = 0
                for j in range(MeshClass.mesh.num_cells()):
                    for k in range(0, len(TracersClass.tracers_in_cells[j])):
                        tracer_no = TracersClass.tracers_in_cells[j][k]
                        if (TracersClass.tracers[tracer_no][8][i] == 1.0): 
                            n_current += 1.0
                
                MPI.barrier(comm)

                n_current = MPI.sum(comm, n_current)
                n_tracers_current.append(n_current)
                n_tracers_orig.append(n_current)

            for i in range(len(materials)):
                composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i], n_tracers_orig, n_tracers_current)
            
            # - Number of tracers -
            tracer_count_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells,  EqClass.number_of_tracers)

            # - Melt fraction -
            if (internal_melting == True):
                melt_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 0, EqClass.xm)
                melt_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 0, EqClass.xm_k)
            
            # - Plastic strain -
            if (plasticity == True):
                scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)

            # - Stress tensor -
            if (elasticity == True):
                stress_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, EqClass.stress_dev_tensor)
                EqClass.stress_dev_inv.assign(project(tensor_2nd_invariant(EqClass.stress_dev_tensor), EqClass.sDG0))
                EqClass.stress_dev_inv_k.assign(project(tensor_2nd_invariant(EqClass.stress_dev_tensor), EqClass.sDG0))

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

        if (use_tracers == True):
            if (initial_topography == True):
                TracersClass.introduce_tracers_unstructured()
            TracersClass.tracer_count_interpolation()
            scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)

            # --- Interpolate composition ---
            n_tracers_orig = []
            n_tracers_current = []
            for i in range(len(materials)):
                n_orig = 0
                for j in range(MeshClass.mesh.num_cells()):
                    for k in range(0, len(TracersClass.tracers_in_cells[j])):
                        tracer_no = TracersClass.tracers_in_cells[j][k]
                        if (TracersClass.tracers[tracer_no][8][i] == 1.0): 
                            n_orig += 1.0
                
                MPI.barrier(comm)

                n_orig = MPI.sum(comm, n_orig)
                n_tracers_orig.append(n_orig)
                n_tracers_current.append(n_orig)

            for i in range(len(materials)):
                composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i], n_tracers_orig, n_tracers_current)

        # Define the temperature even though it will not be solved
        if rank == 0: print("Prescribing initial condition for temperature.")
        EqClass.Temp.assign(project(Expression("Tb - x[1]/h*(Tb-Ts)", Tb=BC_heat_transfer[1][1], Ts=BC_heat_transfer[0][1], h=height, degree=1), ElemClass.sCG2))

        if (solve_energy_problem == True):             
            
            if (init_cond_profile == True):
                EqClass.solve_initial_heat_equation()
            
            if (cos_perturbation == True):
                EqClass.thermal_perturbation()

    
    # --- Termination condition ---
    # - Must not be earlier, need time after reloading -
    if (termination_condition[0] != "initial_condition"):
        if (termination_condition[0] == "time"):
            t_end = termination_condition[1]
        elif (termination_condition[0] == "step"):
            step_end = termination_condition[1]
        elif (termination_condition[0] == "time_and_step"):
            t_end = termination_condition[1]
            step_end = termination_condition[2]

        if (output_frequency[0] == "time"):
            if (reload == False or (reload == True and restart_time == True)):
                time_output  = np.linspace(0, t_end, int(t_end/output_frequency[1] + 1))
            else:
                time_output  = np.linspace(float(t), t_end, int((t_end - float(t))/output_frequency[1] + 1))
        else:
            time_output  = None

    # --- Save initial condition ---
    FilesClass.save_paraview_ini()

    if (termination_condition[0] == "initial_condition"):
        if rank == 0: print("Exiting at initial condition.", flush = True)
        exit()
    
    oscillations_detected = False

    # --- Main time loop ---
    while ((termination_condition[0] == "time" and float(t) < t_end)
        or (termination_condition[0] == "step" and step < step_end)
        or (termination_condition[0] == "time_and_step" and (step < step_end and float(t) < t_end))):
        
        # --- Check cache ---
        if (monitor_cache == True and  rank == 0):
            print("\nCache control", flush = True)
            os.system("ls /home/kihoulou/.cache/dijitso/lib/ | wc -l\n", flush = True)

        # --- Solve Stokes problem ---
        EqClass.solve_Stokes_problem(FilesClass, nonlinear_rheology)

        # --- Update time and step ---
        t.assign(float(t + EqClass.dt))
        step += 1

        MPI.barrier(comm)
        if (rank == 0):
            print("\n----------------------------------------------", flush = True)
            print("\tStep:     ", '{:d}'.format(step), flush = True)
            print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string, flush = True)
            print("\tTime step:", '{:.3e}'.format(float(EqClass.dt/time_units)), time_units_string, flush = True)
            print("----------------------------------------------\n", flush = True)

        # --- Solve heat transfer equation ---
        if (solve_energy_problem == True):
            EqClass.solve_heat_equation()

        # --- Solve motion of the boundaries ---
        if (MeshClass.moving_mesh == True):
            EqClass.solve_topography_evolution(oscillations_detected)

        # --- Advect tracers and move mesh ---
        MPI.barrier(comm)
        # - Advection part I -
        if (use_tracers == True):
            TracersClass.advect_tracers(EqClass.v_k, EqClass.v_mesh, EqClass.dt)
            
        # - Mesh displacement -
        if (MeshClass.moving_mesh == True):
            MeshClass.move_mesh(ElemClass.u_mesh)

        # - Advection part II -
        if (use_tracers == True):
            TracersClass.find_tracers(EqClass.v_k, EqClass.dt)        
            TracersClass.delete_and_find()

            if (TracersClass.only_melt_tracers == False and TracersClass.only_tracked_tracers == False):
                TracersClass.add_tracers(t, step, time_units)
                TracersClass.tracer_count_interpolation()

            # --- Interpolate tracer-carried functions ---
            # - Composition -
            n_tracers_current = []
            for i in range(len(materials)):
                n_current = 0
                for j in range(MeshClass.mesh.num_cells()):
                    for k in range(0, len(TracersClass.tracers_in_cells[j])):
                        tracer_no = TracersClass.tracers_in_cells[j][k]
                        if (TracersClass.tracers[tracer_no][8][i] == 1.0): 
                            n_current += 1.0
                
                MPI.barrier(comm)

                n_current = MPI.sum(comm, n_current)
                n_tracers_current.append(n_current)

            for i in range(len(materials)):
                composition_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, i, EqClass.composition[i], n_tracers_orig, n_tracers_current)
            
            # - Number of tracers -
            tracer_count_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells,  EqClass.number_of_tracers)

            # - Plastic strain -
            scalar_interpolation(MeshClass.mesh, TracersClass.tracers_in_cells, TracersClass.tracers, 6, "ARITM", EqClass.plastic_strain)

        # --- Update stress ---
        if (elasticity == True):
            EqClass.update_stress()
            EqClass.rotate_and_interpolate_stress()

        # --- Postprocessing of numbers ---
        EqClass.top_length.assign(assemble(EqClass.unit_scalar*MeshClass.ds(1)))
        
        # --- Check whether to save results? ---> If yes, save them.
        if (len(save_tracer_trajectory) > 0):
            TracersClass.save_trajectory(t, time_units, step)
        
        step_output, output_now = Output_Timing(step, step_output, t, time_output)
        if (output_now == True):

            save_topography(MeshClass.bound_mesh, name, EqClass.q_ice, EqClass.q_water, step_output)
            if (use_tracers == True):
                # TracersClass.tracer_count_interpolation()
                TracersClass.save_tracers(step_output, EqClass.p_k, EqClass.Temp)

            update_functions_for_output(EqClass, ElemClass)

            TracersClass.rank_interpolation()
            FilesClass.save_paraview(t, time_units)
            FilesClass.save_HDF5(step_output, step, EqClass.dt, t, time_units)

        code_now        = time.time()
        total_time      = (code_now - code_start)/3600.0
        timestep_time   = (code_now - code_now_k)
        code_now_k      = time.time()

        FilesClass.write_statistic(t, step, stat_output, time_units,\
                                q_cond_top  = EqClass.q_cond_top,\
                                q_top       = EqClass.q_top,\
                                v           = EqClass.v_k,\
                                avg_h_bot   = EqClass.h_bot_aver,\
                                h_top_max   = EqClass.h_top,\
                                time        = total_time,\
                                timestep    = timestep_time)

    FilesClass.save_paraview(t, time_units)
    FilesClass.save_HDF5(step_output, step, EqClass.dt, t, time_units)

    if (rank == 0):
        print("\n----------------------------------------------")
        print("\tTermination condition satisfied.")
        print("\tStep:     ", '{:d}'.format(step))
        print("\tTime:     ", '{:.3e}'.format(float(t/time_units)), time_units_string)
        print("----------------------------------------------\n")

def update_functions_for_output(EqClass, ElemClass):
     # --- Postprocessing of functions ---
    if ("density" in paraview_output):
        EqClass.density.assign(project(rho(EqClass.Temp, EqClass.composition, EqClass.xm), ElemClass.sDG0))

    if ("cohesion" in paraview_output):
        EqClass.cohesion.assign(project(cohesion(EqClass.plastic_strain), ElemClass.sDG0))

    if ("mechanisms" in paraview_output):
        ElemClass.mechanisms.assign(project(get_mechanisms(EqClass.Temp, EqClass.v_k, EqClass.stress_dev_inv), ElemClass.sDG0))

    if ("tidal_heating" in paraview_output):
        ElemClass.heating.assign(project(tidal_heating(EqClass.p_k, EqClass.Temp, EqClass.v_k,\
                                EqClass.stress_dev_inv, EqClass.xm, EqClass.composition, EqClass.plastic_strain,\
                                EqClass.step, EqClass.Picard_iter, EqClass.stress_dev_inv_k, EqClass.dt, EqClass.sr_min), ElemClass.sDG0))

    if ("melting_rate" in paraview_output):
        EqClass.melting_rate.assign(project(rho_m*(EqClass.xm - EqClass.xm_k)/EqClass.dt*rho_s/rho_m, ElemClass.sDG0))

    if ("stress_dev_inv" in paraview_output):
        # --- Calculate here only if elasticity and plasticity are off 
        #     otherwise it is already computed in m_equations.py ---
        if (plasticity == False and elasticity == False):
            EqClass.stress_dev_inv.assign(project(2*EqClass.visc*strain_rate_II(EqClass.v_k), ElemClass.sDG0))
    
    if ("strain_rate_inv" in paraview_output):
        if (plasticity == False and elasticity == False):
            EqClass.strain_rate_inv.assign(project(strain_rate_II(EqClass.v_k), EqClass.sDG0))

    if ("viscosity" in paraview_output):
        if (plasticity == False and elasticity == False):
            EqClass.visc.assign(project(eta_ductile(EqClass.Temp, EqClass.v_k, None, EqClass.xm, EqClass.composition, EqClass.step, EqClass.Picard_iter, "mesh"), EqClass.sDG0))
        else:
            EqClass.update_viscosity()

if __name__ == "__main__":
    run_code()