from m_parameters import *
from m_constants import *
import os

def Check_Input_Parameters():

    print("Performing check of input parameters...")

    exit_code = False

    input_par = [1.0, yr, kyr, Myr]
    if (time_units not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid scale_time parameter.")
            print("Possible inputs: 1.0, yr, kyr, Myr")
            print(" ")

    if ((reloading_HDF5 == True or reloading_tracers == True) and reload_name == name):
        exit_code = True
        if (rank == 0):
            print("Reloading file = new file. Cannot overwrite data the code is reading from.")
            print(" ")

    input_par = ["steps", "time"]
    if (output_type not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid 'output_time' parameter:", output_type)
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
    
    if (BC_vel_top not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Top boundary condition for Stokes problem undefined:", BC_vel_top)
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_vel_bot not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Bottom boundary condition for Stokes problem undefined.", BC_vel_bot)
            print("Possible inputs:", bc_list)
            print(" ")

    bc_list = ["free_slip", "no_slip", "velocity", "velocity_x", "velocity_y"]

    if (BC_vel_left not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Left boundary condition for Stokes problem undefined.", BC_vel_left)
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_vel_right not in bc_list):
        exit_code = True
        if (rank == 0):
            print("Right boundary condition for Stokes problem undefined.", BC_vel_right)
            print("Possible inputs:", bc_list)
            print(" ")

    if (BC_vel_left == "free_surface" or BC_vel_right == "free_surface"):
        exit_code = True
        if (rank == 0):
            print("No free surface boudnary condition for side boundaries.")
            print("Possible inputs:", bc_list)
            print(" ")

    bc_temp_top_list = ["temperature", "heat_flux", "radiation"]
    
    if (BC_T_top not in bc_temp_top_list):
        exit_code = True
        if (rank == 0):
            print("Top boundary condition for Energy balance undefined:", BC_T_top)
            print("Possible inputs:", bc_temp_top_list)
            print(" ")
    
    bc_temp_list = ["temperature", "heat_flux"]
    if (BC_T_bot not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Bottom boundary condition for Energy balance undefined:", BC_T_bot)
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_T_left not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Left boundary condition for Energy balance undefined:", BC_T_left)
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if (BC_T_right not in bc_temp_list):
        exit_code = True
        if (rank == 0):
            print("Right boundary condition for Energy balance undefined:", BC_T_right)
            print("Possible inputs:", bc_temp_list)
            print(" ")

    if ((BC_T_top == "temperature" and T_top < 0.0)\
        or (BC_T_bot == "temperature" and T_bot < 0.0)\
        or (BC_T_right == "temperature" and T_right < 0.0)\
        or (BC_T_left == "temperature" and T_left < 0.0)):

        exit_code = True
        if (rank == 0):
            print("Temperature boundary condition cannot be less than 0 K.")
            print(" ")

    if ((BC_T_top == "heat_flux" or BC_T_top == "radiation") and BC_T_bot == "heat_flux" and BC_T_left == "heat_flux" and BC_T_right == "heat_flux"):

        exit_code = True
        if (rank == 0):
            print("At least one Dirichlet boundary condition needs to be prescribed for the Energy balance.")
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
    if (timestep_strategy not in input_par):
        exit_code = True
        if(rank == 0):
            print("Invalid time step strategy:", timestep_strategy)
            print("Possible inputs:", input_par)
            print(" ")

    if (exit_code == True):
        exit()
    else:
        print("Check done.\n")

def SaveSourceCode():
    if (rank == 0): 
        os.system("cp  main*.py data_" + name + "/source_code")
        os.system("cp  m_*.py data_" + name + "/source_code")	