from m_parameters import *
from m_constants import *

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKCYAN = '\033[96m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

def check_input_parameters():

    if (rank == 0):
        print("Performing check of input parameters...")

    exit_code = False
    count = 0

    if ("mechanisms" in paraview_output and (viscosity_type == "temp-dep" or viscosity_type == "constant")):
        exit_code = True
        if (rank == 0):
            count += 1
            print("'Mechanisms' only defined for stress-dependent rheology.")
            print(" ")

    if (len(save_tracer_trajectory) > 0):
        for i in range(len(save_tracer_trajectory)):
            if (isinstance(save_tracer_trajectory[i][2], str) == False):
                exit_code = True
                if (rank == 0):
                    count += 1
                    print("'save_tracer_trajectory' at index ", i ,"requires a string ID.")
                    print(" ")

    input_par = ["-", "s", "hr", "day", "yr", "kyr", "Myr", "Gyr"]
    if (time_units_string not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid scale_time parameter.")
            print("Possible inputs: '-', 's', 'hr', 'day', 'yr', 'kyr', 'Myr', 'Gyr'")
            print(" ")

    if ((reload == True) and reload_name == name):
        exit_code = True
        if (rank == 0):
            count += 1
            print(f"{FAIL}Reloading file = new file.{ENDC}" + " Cannot overwrite data the code is reading from.")
            print(" ")
    
    input_par = ["step", "time", "time_and_step," "initial_condition"]
    if (termination_condition[0] not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid 'termination_condition' parameter:", termination_condition[0])
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["steps", "time"]
    if (output_frequency[0] not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid 'output_frequency' parameter:", output_frequency[0])
            print("Possible inputs:", input_par)
            print(" ")

    if (reload == False):
        if (len(refinement)==0 or len(refinement)%4==0):
            pass
        else:
            exit_code = True
            if(rank == 0):
                count += 1
                print("Invalid mesh refinement input:", refinement)
                print("Input should follow this format: [x_left, x_right, y_bottom, y_top, etc.]")
                print(" ")

    bc_list = ["free_slip", "no_slip", "free_surface", "velocity", "velocity_x", "velocity_y"]
    
    if (BC_Stokes_problem[0][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Top boundary condition for Stokes problem undefined:", BC_Stokes_problem[0][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[1][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Bottom boundary condition for Stokes problem undefined.", BC_Stokes_problem[1][0])
            print("Possible inputs:", bc_list)
            print(" ")

    bc_list = ["free_slip", "no_slip", "velocity", "velocity_x", "velocity_y"]

    if (BC_Stokes_problem[2][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Left boundary condition for Stokes problem undefined.", BC_Stokes_problem[2][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[3][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Right boundary condition for Stokes problem undefined.", BC_Stokes_problem[3][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[2][0] == "free_surface" or BC_Stokes_problem[3][0] == "free_surface"):
        exit_code = True
        if (rank == 0):
            count += 1
            print("No free surface boudnary condition for side boundaries.")
            print("Possible inputs:", bc_list)
            print(" ")

    bc_temp_top_list = ["temp", "heat_flux", "radiation"]
    
    if (BC_heat_transfer[0][0] not in bc_temp_top_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print(f"Top boundary condition for Energy balance undefined:{WARNING}", BC_heat_transfer[0][0], f"{ENDC}")
            print(f"Possible inputs:{OKGREEN}", bc_temp_top_list, f"{ENDC}")
            print(" ")
    
    bc_temp_list = ["temp", "heat_flux"]
    if (BC_heat_transfer[1][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Bottom boundary condition for Energy balance undefined:", BC_heat_transfer[1][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_heat_transfer[2][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Left boundary condition for Energy balance undefined:", BC_heat_transfer[2][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_heat_transfer[3][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            count += 1
            print("Right boundary condition for Energy balance undefined:", BC_heat_transfer[3][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if ((BC_heat_transfer[0][0] == "temp"    and BC_heat_transfer[0][1] < 0.0)\
        or (BC_heat_transfer[1][0] == "temp" and BC_heat_transfer[1][1] < 0.0)\
        or (BC_heat_transfer[2][0] == "temp" and BC_heat_transfer[2][1] < 0.0)\
        or (BC_heat_transfer[3][0] == "temp" and BC_heat_transfer[3][1] < 0.0)):

        exit_code = True
        if (rank == 0):
            count += 1
            print("Temperature boundary condition cannot be less than 0 K.")
            print(" ")

    if (BC_heat_transfer[0][0] != "temp" and BC_heat_transfer[1][0] != "temp" and BC_heat_transfer[2][0] != "temp" and BC_heat_transfer[3][0] != "temp"):

        exit_code = True
        if (rank == 0):
            count += 1
            print("At least one Dirichlet boundary condition needs to be prescribed for the heat transfer problem.")
            print(" ")

    input_par = ["constant", "temp-dep", "GK_2001", "composition"]
    if (viscosity_type not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print(f"Invalid viscosity type:{WARNING}", viscosity_type, f"{ENDC}")
            print(f"Possible inputs:{OKGREEN}", input_par, f"{ENDC}")
            print(" ")

    input_par = ["Euler", "RK2", "RK4"]
    if (integration_method not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid integration method type:", integration_method)
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["Maxwell", "Andrade", "none"]
    if (heating_model not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid tidal heating model:", heating_model)
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["domain", "cell", "constant", "convective"]
    if (time_step_strategy not in input_par):
        exit_code = True
        if(rank == 0):
            count += 1
            print("Invalid time step strategy:", time_step_strategy)
            print("Possible inputs:", input_par)
            print(" ")
    
    # --- Ensuring that the parameters needed for restarting are always written ---
    # ---- Field functions ---
    if (solve_energy_problem == True and "temperature" not in paraview_output):
        paraview_output.append("temperature")
        if (rank == 0):
            print("Adding 'temperature' into paraview_output list.")

    if (BC_Stokes_problem[0][0] == "free_surface" and "topography_top" not in paraview_output):
        paraview_output.append("topography_top")
        if (rank == 0):
            print("Adding 'topography_top' into paraview_output list.")

    if (BC_Stokes_problem[1][0] == "free_surface" and "topography_bottom" not in paraview_output):
        paraview_output.append("topography_bottom")
        if (rank == 0):
            print("Adding 'topography_bottom' into paraview_output list.")

    if (internal_melting == True and "melt_fraction" not in paraview_output):
        paraview_output.append("melt_fraction")
        if (rank == 0):
            print("Adding 'melt_fraction' into paraview_output list.")

    if (plasticity == True and "pressure" not in paraview_output):
        paraview_output.append("pressure")
        if (rank == 0):
            print("Adding 'pressure' into paraview_output list.")

    # ---- Tracers ---
    # --- The headers need to be embedded so that they can be recognized when reloading ---
    # --- Loop over the ones that are in tracers_output --- 
    for arg in tracers_output:
        if (arg == "rank"):
            tracers_header.append("rank")

        if (arg == "dev_stress_xx"):
            tracers_header.append("s_xx (Pa)")
        
        if (arg == "dev_stress_xz"):
            tracers_header.append("s_xz (Pa)")

        if (arg == "plastic_strain"):
            tracers_header.append("eps_p (-)")

        if (arg == "original_y"):
            tracers_header.append("y_orig (m)")

        if (arg == "melt_fraction"):
            tracers_header.append("xm (-)")

        if (arg == "origin"):
            tracers_header.append("origin")

        if (arg == "ID"):
            tracers_header.append("ID")

    # --- Now add the mandatory ones, if not there already ---
    if (internal_melting == True and "melt_fraction" not in tracers_output):
        tracers_output.append("melt_fraction")
        tracers_header.append("xm (-)")
        if (rank == 0):
            print("Adding 'melt_fraction' into tracers_output list.")

    if (plasticity == True and "plastic_strain" not in tracers_output):
        tracers_output.append("plastic_strain")
        tracers_header.append("eps_p (-)")
        if (rank == 0):
            print("Adding 'plastic_strain' into tracers_output list.")

    if (plasticity == True):
        tracers_output.append("origin")
        tracers_header.append("origin")
        if (rank == 0):
            print("Adding 'origin' into tracers_output list.")

    if (elasticity == True and "dev_stress_xx" not in tracers_output):
        tracers_output.append("dev_stress_xx")
        tracers_header.append("s_xx (Pa)")
        if (rank == 0):
            print("Adding 'dev_stress_xx' into tracers_output list.")

    if (elasticity == True and "dev_stress_xz" not in tracers_output):
        tracers_output.append("dev_stress_xz")
        tracers_header.append("s_xz (Pa)")
        if (rank == 0):
            print("Adding 'dev_stress_xz' into tracers_output list.")

    # --- This one has to be the last because it can produce multiple columns ---
    if (len(materials) > 0 and "composition" not in tracers_output):
        tracers_output.append("composition")
        tracers_header.append("composition")
        if (rank == 0):
            print("Adding 'composition' into tracers_output list.")

    if (exit_code == True):
        if (rank == 0):
            if (count == 1):
                print(f"{FAIL}{BOLD}There is " + str(count) + f" incorrect input.{ENDC}\n")
            else:
                print(f"{FAIL}{BOLD}There are " + str(count) + f" incorrect inputs.{ENDC}\n")
        exit()
    else:
        if (rank == 0):
            print("No errors in input were found.")