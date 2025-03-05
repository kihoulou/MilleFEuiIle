from dolfin import *
from m_parameters import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

def apply_temperature_BC(sCG2, boundary_parts):
    bc_T_top = DirichletBC(sCG2, T_top, boundary_parts,1)
    bc_T_bot = DirichletBC(sCG2, T_bot, boundary_parts,2)
    bc_T_left = DirichletBC(sCG2, T_left, boundary_parts,3)
    bc_T_right = DirichletBC(sCG2, T_right, boundary_parts,4)

    bc_temp = []

    if (BC_T_top == "temperature"):
        bc_temp.append(bc_T_top)

    if (BC_T_bot == "temperature"):
        bc_temp.append(bc_T_bot)

    if (BC_T_left == "temperature"):
        bc_temp.append(bc_T_left)

    if (BC_T_right == "temperature"):
        bc_temp.append(bc_T_right)
    
    return bc_temp

class Point_Fixed_Pressure(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) and near(x[1], height)

def apply_velocity_BC(V, boundary_parts, top_left):

    bc_stokes = []

    # --- Free slip conditions ---
    bc_free_slip_top = DirichletBC(V.sub(1).sub(1), Constant(0.0), boundary_parts,1)
    bc_free_slip_bot = DirichletBC(V.sub(1).sub(1), Constant(0.0), boundary_parts,2)
    bc_free_slip_left = DirichletBC(V.sub(1).sub(0), Constant(0.0), boundary_parts,3)
    bc_free_slip_right = DirichletBC(V.sub(1).sub(0), Constant(0.0), boundary_parts,4)

    # --- No slip conditions ---
    bc_no_slip_top = DirichletBC(V.sub(1), Constant((0.0,0.0)),boundary_parts,1)
    bc_no_slip_bot = DirichletBC(V.sub(1), Constant((0.0,0.0)),boundary_parts,2)
    bc_no_slip_left = DirichletBC(V.sub(1), Constant((0.0,0.0)),boundary_parts,3)
    bc_no_slip_right = DirichletBC(V.sub(1), Constant((0.0,0.0)),boundary_parts,4)

    # --- General velocity conditions ---
    bc_velocity_top_x = DirichletBC(V.sub(1).sub(0), velocity_top_x, boundary_parts,1)
    bc_velocity_bot_x = DirichletBC(V.sub(1).sub(0), velocity_bot_x, boundary_parts,2)
    bc_velocity_left_x = DirichletBC(V.sub(1).sub(0), velocity_left_x, boundary_parts,3)
    bc_velocity_right_x = DirichletBC(V.sub(1).sub(0), velocity_right_x, boundary_parts,4)

    bc_velocity_top_y = DirichletBC(V.sub(1).sub(1), velocity_top_y, boundary_parts,1)
    bc_velocity_bot_y = DirichletBC(V.sub(1).sub(1), velocity_bot_y, boundary_parts,2)
    bc_velocity_left_y = DirichletBC(V.sub(1).sub(1), velocity_left_y, boundary_parts,3)
    bc_velocity_right_y = DirichletBC(V.sub(1).sub(1), velocity_right_y, boundary_parts,4)

    # --- Fixed pressure condition ---
    bc_pres = DirichletBC(V.sub(0), Constant(0.0), top_left, method = "pointwise") 
    

    # --- Top boundary ---
    if (BC_vel_top == "free_slip"):
        bc_stokes.append(bc_free_slip_top)

    elif (BC_vel_top == "no_slip"):
        bc_stokes.append(bc_no_slip_top)

    elif (BC_vel_top == "free_surface"):
        pass

    elif (BC_vel_top == "velocity_x"):
        bc_stokes.append(bc_velocity_top_x)

    elif (BC_vel_top == "velocity_y"):
        bc_stokes.append(bc_velocity_top_y)

    elif (BC_vel_top == "velocity"):
        bc_stokes.append(bc_velocity_top_x)
        bc_stokes.append(bc_velocity_top_y)

    else:
        if (rank == 0):
            print("Top boundary condition undefined. Exiting.")
        exit()

    # --- Bottom boundary ---
    if (BC_vel_bot == "free_slip"):
        bc_stokes.append(bc_free_slip_bot)

    elif (BC_vel_bot == "no_slip"):
        bc_stokes.append(bc_no_slip_bot)

    elif (BC_vel_bot == "free_surface"):
        pass

    elif (BC_vel_bot == "velocity_x"):
        bc_stokes.append(bc_velocity_bot_x)

    elif (BC_vel_bot == "velocity_y"):
        bc_stokes.append(bc_velocity_bot_y)

    elif (BC_vel_bot == "velocity"):
        bc_stokes.append(bc_velocity_bot_x)
        bc_stokes.append(bc_velocity_bot_y)

    else:
        if (rank == 0):
            print("Bottom boundary condition undefined. Exiting.")
        exit()
    
    # --- Left boundary ---
    if (BC_vel_left == "free_slip"):
        bc_stokes.append(bc_free_slip_left)

    elif (BC_vel_left == "no_slip"):
        bc_stokes.append(bc_no_slip_left)

    elif (BC_vel_left == "free_surface"):
        pass

    elif (BC_vel_left == "velocity_x"):
        bc_stokes.append(bc_velocity_left_x)

    elif (BC_vel_left == "velocity_y"):
        bc_stokes.append(bc_velocity_left_y)

    elif (BC_vel_left == "velocity"):
        bc_stokes.append(bc_velocity_left_x)
        bc_stokes.append(bc_velocity_left_y)

    else:
        if (rank == 0):
            print("Left boundary condition undefined. Exiting.")
        exit()

    # --- Right boundary ---
    if (BC_vel_right == "free_slip"):
        bc_stokes.append(bc_free_slip_right)

    elif (BC_vel_right == "no_slip"):
        bc_stokes.append(bc_no_slip_right)

    elif (BC_vel_right == "free_surface"):
        pass

    elif (BC_vel_right == "velocity_x"):
        bc_stokes.append(bc_velocity_right_x)

    elif (BC_vel_right == "velocity_y"):
        bc_stokes.append(bc_velocity_right_y)

    elif (BC_vel_right == "velocity"):
        print("appluing ")
        bc_stokes.append(bc_velocity_right_x)
        bc_stokes.append(bc_velocity_right_y)

    else:
        if (rank == 0):
            print("Right boundary condition undefined. Exiting.")
        exit()
    
    if (BC_vel_top != "free_surface" and BC_vel_bot != "free_surface" and BC_vel_left != "free_surface" and BC_vel_right != "free_surface"):
        bc_stokes.append(bc_pres)

    return bc_stokes