from m_parameters import *
from m_constants import *
import os

def check_input_parameters():
    """Performs a check on all string-type inputs in ``m_parameters`` module before running the code.

    :returns: Parameters that have an incorrect input (if there are any) and their admissible values.

    .. note:: Only the string-type parameters are being checked.
    """

    print("Performing check of input parameters...")

    exit_code = False

    input_par = [1.0, yr, kyr, Myr]
    if (time_units not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid scale_time parameter.")
            print("Possible inputs: 1.0, yr, kyr, Myr")
            print(" ")

    if ((reload_HDF5 == True or reload_tracers == True) and reload_name == name):
        exit_code = True
        if (rank == 0):
            print("Reloading file = new file. Cannot overwrite data the code is reading from.")
            print(" ")

    input_par = ["steps", "time"]
    if (output_frequency[0] not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid 'output_time' parameter:", output_frequency[0])
            print("Possible inputs:", input_par)
            print(" ")

    if (loading_mesh == False):
        if (len(refinement)==0 or len(refinement)%4==0):
            pass
        else:
            exit_code = True
            if(rank == 0):
                print("Invalid mesh refinement input:", refinement)
                print("Input should follow this format: [x_left, x_right, y_bottom, y_top, etc.]")
                print(" ")

    bc_list = ["free_slip", "no_slip", "free_surface", "velocity", "velocity_x", "velocity_y"]
    
    if (BC_Stokes_problem[0][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Top boundary condition for Stokes problem undefined:", BC_Stokes_problem[0][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[1][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Bottom boundary condition for Stokes problem undefined.", BC_Stokes_problem[1][0])
            print("Possible inputs:", bc_list)
            print(" ")

    bc_list = ["free_slip", "no_slip", "velocity", "velocity_x", "velocity_y"]

    if (BC_Stokes_problem[2][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Left boundary condition for Stokes problem undefined.", BC_Stokes_problem[2][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[3][0] not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Right boundary condition for Stokes problem undefined.", BC_Stokes_problem[3][0])
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_Stokes_problem[2][0] == "free_surface" or BC_Stokes_problem[3][0] == "free_surface"):
        exit_code = True
        if (rank == 0):
            print("No free surface boudnary condition for side boundaries.")
            print("Possible inputs:", bc_list)
            print(" ")

    bc_temp_top_list = ["temp", "heat_flux", "radiation"]
    
    if (BC_heat_transfer[0][0] not in bc_temp_top_list):
        exit_code = True
        if (rank == 0):
            print("Top boundary condition for Energy balance undefined:", BC_heat_transfer[0][0])
            print("Possible inputs:", bc_temp_top_list)
            print(" ")
    
    bc_temp_list = ["temp", "heat_flux"]
    if (BC_heat_transfer[1][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Bottom boundary condition for Energy balance undefined:", BC_heat_transfer[1][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_heat_transfer[2][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Left boundary condition for Energy balance undefined:", BC_heat_transfer[2][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_heat_transfer[3][0] not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Right boundary condition for Energy balance undefined:", BC_heat_transfer[3][0])
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if ((BC_heat_transfer[0][0] == "temp"    and BC_heat_transfer[0][1] < 0.0)\
        or (BC_heat_transfer[1][0] == "temp" and BC_heat_transfer[1][1] < 0.0)\
        or (BC_heat_transfer[2][0] == "temp" and BC_heat_transfer[2][1] < 0.0)\
        or (BC_heat_transfer[3][0] == "temp" and BC_heat_transfer[3][1] < 0.0)):

        exit_code = True
        if (rank == 0):
            print("Temperature boundary condition cannot be less than 0 K.")
            print(" ")

    if (BC_heat_transfer[0][0] != "temp" and BC_heat_transfer[1][0] != "temp" and BC_heat_transfer[2][0] != "temp" and BC_heat_transfer[3][0] != "temp"):

        exit_code = True
        if (rank == 0):
            print("At least one Dirichlet boundary condition needs to be prescribed for the heat transfer problem.")
            print(" ")

    input_par = ["constant", "temp-dep", "GK_2001", "composition"]
    if (viscosity_type not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid viscosity type:", viscosity_type)
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["Euler", "RK2", "RK4"]
    if (integration_method not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid integration method type:", integration_method)
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["Maxwell", "Andrade", "none"]
    if (heating_model not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid tidal heating model:", heating_model)
            print("Possible inputs:", input_par)
            print(" ")

    input_par = ["domain", "cell", "constant"]
    if (time_step_strategy not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid time step strategy:", time_step_strategy)
            print("Possible inputs:", input_par)
            print(" ")

    if (exit_code == True):
        exit()
    else:
        print("Check done.\n")