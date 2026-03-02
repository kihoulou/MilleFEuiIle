from dolfin import *
from m_parameters_docs import *
from m_constants import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

def apply_temperature_BC(sCG2, boundary_parts):

    bc_temp = []

    # --- Top boundary ---
    if (BC_heat_transfer[0][0] == "temp"):
        bc_temp.append(DirichletBC(sCG2, BC_heat_transfer[0][1], boundary_parts, 1))

    # --- Bottom boundary ---
    if (BC_heat_transfer[1][0] == "temp"):
        bc_temp.append(DirichletBC(sCG2, BC_heat_transfer[1][1], boundary_parts, 2))

    # --- Left boundary ---
    if (BC_heat_transfer[2][0] == "temp"):
        bc_temp.append(DirichletBC(sCG2, BC_heat_transfer[2][1], boundary_parts, 3))

    # --- Right boundary ---
    if (BC_heat_transfer[3][0] == "temp"):
        bc_temp.append(DirichletBC(sCG2, BC_heat_transfer[3][1], boundary_parts, 4))
    
    return bc_temp

class Point_Fixed_Pressure(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], length) and near(x[1], height)

def apply_velocity_BC(V, boundary_parts, top_left, t):

    bc_stokes = []
    
    # --- Top boundary ---
    if (BC_Stokes_problem[0][0] == "free_slip"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(0.0), boundary_parts, 1))

    elif (BC_Stokes_problem[0][0] == "no_slip"):
        bc_stokes.append(DirichletBC(V.sub(1), Constant((0.0, 0.0)), boundary_parts, 1))

    elif (BC_Stokes_problem[0][0] == "free_surface"):
        pass

    elif (BC_Stokes_problem[0][0] == "velocity_x"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[0][1]), boundary_parts, 1))

    elif (BC_Stokes_problem[0][0] == "velocity_y"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[0][1]), boundary_parts, 1))

    elif (BC_Stokes_problem[0][0] == "velocity"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[0][1]), boundary_parts, 1))
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[0][2]), boundary_parts, 1))

    else:
        if (rank == 0):
            print("Top boundary condition undefined. Exiting.")
        exit()

    def left_x_time_dependent(t):
        # if (float(t) <= 2*kyr):
        #     return -10e3/Myr
        # else:
        #     return 0.0
        return -10e3/Myr*sin(2*pi*t/(6*Myr))
                
    def right_x_time_dependent(t):
        return +10e3/Myr*sin(2*pi*t/(6*Myr))
    
    def bottom_y_time_dependent(t):
        return length/(2*height)*(-left_x_time_dependent(t) + right_x_time_dependent(t))
    

    # --- Bottom boundary ---
    if (BC_Stokes_problem[1][0] == "free_slip"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(0.0), boundary_parts, 2))

    elif (BC_Stokes_problem[1][0] == "no_slip"):
        bc_stokes.append(DirichletBC(V.sub(1), Constant((0.0,0.0)), boundary_parts, 2))

    elif (BC_Stokes_problem[1][0] == "free_surface"):
        pass

    elif (BC_Stokes_problem[1][0] == "velocity_x"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[1][1]), boundary_parts, 2))

    elif (BC_Stokes_problem[1][0] == "velocity_y"):
        
        if (BC_Stokes_problem[2][1] == "time_dependent"):
            bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant((bottom_y_time_dependent(t))), boundary_parts, 2))
        else:
            bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[1][1]), boundary_parts, 2))

    elif (BC_Stokes_problem[1][0] == "velocity"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[1][1]), boundary_parts, 2))
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[1][2]), boundary_parts, 2))

    else:
        if (rank == 0):
            print("Bottom boundary condition undefined. Exiting.")
        exit()
    
    # --- Left boundary ---
    if (BC_Stokes_problem[2][0] == "free_slip"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(0.0), boundary_parts, 3))

    elif (BC_Stokes_problem[2][0] == "no_slip"):
        bc_stokes.append(DirichletBC(V.sub(1), Constant((0.0,0.0)), boundary_parts, 3))

    elif (BC_Stokes_problem[2][0] == "free_surface"):
        pass

    elif (BC_Stokes_problem[2][0] == "velocity_x"):
        if (BC_Stokes_problem[2][1] == "time_dependent"):
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant((left_x_time_dependent(t))), boundary_parts, 3))
        else:
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[2][1]), boundary_parts, 3))

    elif (BC_Stokes_problem[2][0] == "velocity_y"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[2][1]), boundary_parts, 3))

    elif (BC_Stokes_problem[2][0] == "velocity"):
        if (BC_Stokes_problem[2][1] == "time_dependent"):
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant((left_x_time_dependent(t))), boundary_parts, 3))
        else:
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[2][1]), boundary_parts, 3))
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[2][2]), boundary_parts, 3))

    else:
        if (rank == 0):
            print("Left boundary condition undefined. Exiting.")
        exit()

    # --- Right boundary ---
    if (BC_Stokes_problem[3][0] == "free_slip"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(0.0), boundary_parts, 4))

    elif (BC_Stokes_problem[3][0] == "no_slip"):
        bc_stokes.append(DirichletBC(V.sub(1), Constant((0.0,0.0)), boundary_parts, 4))

    elif (BC_Stokes_problem[3][0] == "free_surface"):
        pass

    elif (BC_Stokes_problem[3][0] == "velocity_x"):
        if (BC_Stokes_problem[3][1] == "time_dependent"):
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant((right_x_time_dependent(t))), boundary_parts, 4))
        else:
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[3][1]), boundary_parts, 4))

    elif (BC_Stokes_problem[3][0] == "velocity_y"):
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[3][1]), boundary_parts, 4))

    elif (BC_Stokes_problem[3][0] == "velocity"):
        if (BC_Stokes_problem[3][1] == "time_dependent"):
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant((right_x_time_dependent(t))), boundary_parts, 4))
        else:
            bc_stokes.append(DirichletBC(V.sub(1).sub(0), Constant(BC_Stokes_problem[3][1]), boundary_parts, 4))
        bc_stokes.append(DirichletBC(V.sub(1).sub(1), Constant(BC_Stokes_problem[3][2]), boundary_parts, 4))

    else:
        if (rank == 0):
            print("Right boundary condition undefined. Exiting.")
        exit()
    
    if (BC_Stokes_problem[0][0] != "free_surface" and BC_Stokes_problem[1][0] != "free_surface"):
        bc_stokes.append(DirichletBC(V.sub(0), Constant(0.0), top_left, method = "pointwise"))

    return bc_stokes