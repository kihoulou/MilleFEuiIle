from dolfin import *
from m_parameters import *
import time

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

class Tracers:
    def __init__(self, MeshClass, ElemClass, FilesClass):        
        self.mesh = MeshClass.mesh
        self.moving_mesh = MeshClass.moving_mesh
        
        self.Temp = ElemClass.Temp
        
        self.stress_dev_tensor  = ElemClass.stress_dev_tensor
        self.plastic_strain     = ElemClass.plastic_strain
        self.ocean_material     = ElemClass.ocean_material
        self.surface_material   = ElemClass.surface_material
        self.h_top              = ElemClass.h_top
        self.xm                 = ElemClass.xm
        self.number_of_tracers  = ElemClass.number_of_tracers
        self.composition        = ElemClass.composition
        self.height_fraction    = ElemClass.height_fraction
        self.mesh_ranks    = ElemClass.mesh_ranks

        self.name = FilesClass.name

        # --- Determine whether to use tracers or not ---
        self.use_tracers = False
        if (plasticity == True or elasticity == True or internal_melting == True or len(materials) > 0):
            self.use_tracers = True
            
        # --- Only melt tracers exception ---
        self.only_melt_tracers = False
        if (internal_melting == True):
            if (plasticity == True or elasticity == True or len(materials) > 0):
                self.only_melt_tracers = False
            else:
                self.only_melt_tracers = True
                
        # --- If the tracers are needed, initialize them ---
        if (self.use_tracers == True):
            self.tracers = []
            self.vacancy = []
            self.tracers_in_cells = []

            self.moving_tracers = []
            self.tracers_to_find = []
            
            if (reload_tracers == False):
                if (self.only_melt_tracers == False):
                    self.introduce_tracers()
                else:
                    for j in range(self.mesh.num_cells()):
                        self.tracers_in_cells.append([])
                    
            else:
                self.load_tracers()

    def save_tracers(self, step):
        file = open("data_" + self.name + "/tracers/step_" + str(step) + ".dat","a")
        self.save_header_tracer(file)
        for j in range (0, len(self.tracers)): 
            if (j not in self.vacancy):
                self.save_indiv_tracer(file, j,\
                            rank            = self.tracers[j][3],\
                            dev_stress_xx   = self.tracers[j][4],\
                            dev_stress_xz   = self.tracers[j][5],\
                            plastic_strain  = self.tracers[j][6],\
                            ocean_material  = self.tracers[j][7],\
                            surface_material= self.tracers[j][8],\
                            original_depth  = self.tracers[j][9],\
                            composition     = self.tracers[j][10],\
                            melt_fraction   = self.tracers[j][11][0],\
                            origin          = self.tracers[j][12],\
                            id              = self.tracers[j][13])
        file.close()
    
    def save_header_tracer(self, file):
        file.write((2*"%s\t\t")%("x_pos (m)", "y_pos (m)"))   

        for arg in Tracers_header:
            file.write(("%s\t\t")%(arg))
        file.write("\n")

    def save_indiv_tracer(self, file, j, **kwargs):
        file.write((2*"%.7E\t")%(self.tracers[j][0], self.tracers[j][1]))   

        for arg in Tracers_Output:
            for key in kwargs:
                if (arg == key):
                    if (arg == "composition"):
                        for i in range(len(kwargs[key])):
                            file.write(("%.5E\t")%(kwargs[key][i]))
                    else:
                        file.write(("%.5E\t")%(kwargs[key]))
        file.write("\n")

    def rank_interpolation(self):
        ranks = []
        for j in range(self.mesh.num_cells()):
            ranks.append(rank)

        ranks = np.array(ranks)                                                            
        self.mesh_ranks.vector().set_local(ranks)   

    def tracer_count_interpolation(self):
        ranks = []
        for j in range(self.mesh.num_cells()):
            ranks.append(len(self.tracers_in_cells[j]))

        ranks = np.array(ranks)                                                            
        self.number_of_tracers.vector().set_local(ranks)

    def introduce_tracers(self):
        tracer_no = 0

        for j in range(self.mesh.num_cells()):
            centroid = Cell(self.mesh, j).midpoint()
            cell_j = Cell(self.mesh,j).get_vertex_coordinates() # gives [x_1, y_1, x_2, y_2, x_3, y_3]

            self.tracers_in_cells.append([])

            # --- If the cell lays in a region that should be without tracers, skip the rest ---
            skip_cell = False
            if (empty_cells_allowed == True):
                for empty in empty_cells_region:
                    if (empty[0] == "rectangle"):
                            if ((empty[1] <= centroid.x() <= empty[2]) and (empty[3] <= centroid.y() <= empty[4])):
                                skip_cell = True

            if (skip_cell == True):
                continue

            if ((int(cell_j[0]) != int(cell_j[2])) and (int(cell_j[2]) != int(cell_j[4])) and (int(cell_j[0]) != int(cell_j[4]))):
                order = [int(cell_j[1]), int(cell_j[3]), int(cell_j[5])]
                order.sort()
                if (order[0] == order[1]):
                    orientation = "UP"
                if (order[1] == order[2]):
                    orientation = "DOWN"

            if ((int(cell_j[1]) != int(cell_j[3])) and (int(cell_j[3]) != int(cell_j[5])) and (int(cell_j[1]) != int(cell_j[5]))):
                order = [int(cell_j[0]),int(cell_j[2]),int(cell_j[4])]
                order.sort()
                if (order[0] == order[1]):
                    orientation = "RIGHT"
                if (order[1] == order[2]):
                    orientation = "LEFT"

            if (((int(cell_j[0]) == int(cell_j[2])) or (int(cell_j[0]) == int(cell_j[4])) or (int(cell_j[2]) == int(cell_j[4])))\
                and ((int(cell_j[1]) == int(cell_j[3])) or (int(cell_j[1]) == int(cell_j[5])) or (int(cell_j[3]) == int(cell_j[5])))):
                order = [int(cell_j[0]),int(cell_j[2]),int(cell_j[4])]
                order.sort()

                if (order[0] == order[1]):
                    orientation = "DIAG_LEFT"
                if (order[1] == order[2]):
                    orientation = "DIAG_RIGHT"

            if (orientation == "DOWN" or orientation == "RIGHT" or orientation == "DIAG_RIGHT" or orientation == "DIAG_LEFT"):
                origin_x = min(cell_j[0],min(cell_j[2],cell_j[4]))
                origin_y = max(cell_j[1],max(cell_j[3],cell_j[5]))
                sign = 1.0

            if (orientation == "UP" or orientation == "LEFT"):
                origin_x = max(cell_j[0],max(cell_j[2],cell_j[4]))
                origin_y = min(cell_j[1],min(cell_j[3],cell_j[5]))
                sign = -1.0

            if (orientation == "UP" or orientation == "DOWN" or orientation == "DIAG_RIGHT" or orientation == "DIAG_LEFT"):
                x_extent = max(cell_j[0],max(cell_j[2],cell_j[4])) - min(cell_j[0],min(cell_j[2],cell_j[4]))
                z_extent = max(cell_j[1],max(cell_j[3],cell_j[5])) - min(cell_j[1],min(cell_j[3],cell_j[5]))

            if (orientation == "LEFT" or orientation == "RIGHT"):
                x_extent = max(cell_j[0],max(cell_j[2],cell_j[4])) - min(cell_j[0],min(cell_j[2],cell_j[4]))
                z_extent = max(cell_j[1],max(cell_j[3],cell_j[5])) - min(cell_j[1],min(cell_j[3],cell_j[5]))
                
            lim = int(sqrt(tracers_per_cell*4.0))

            # self.tracers_in_cells.append([])

            dx = x_extent/(lim + 1.0)
            dz = z_extent/(lim + 1.0)

            for k in range(0, lim + 1):
                for l in range(0, lim + 1):

                    x_pos = origin_x + sign*k*dx + dx/2.0
                    y_pos = origin_y - sign*l*dz - dz/2.0

                    if (Cell(self.mesh, j).contains(Point(x_pos, y_pos)) == True):

                        # --- Assign a layer of near-surface ice ---
                        if (y_pos > height - sm_thickness):
                            surf_mat = 1.0
                        else:
                            surf_mat = 0.0

                        # --- Assign a layer of freshly crystallized ice ---
                        if (y_pos < om_thickness):
                            ocean_mat = 1.0
                        else:
                            ocean_mat = 0.0

                        # --- Assign chemical composition ---
                        if (len(materials) == 0):
                            comp = []

                        else:
                            comp = [0.0]*len(materials)
                            i = 0

                            for mat in materials:
                                if (mat[0] == "circle"):
                                    if (sqrt((x_pos - mat[1])**2 + (y_pos - mat[2])**2) < mat[3]):
                                        for ii in range(0, i): # Erase previous composition
                                            comp[ii] = 0.0
                                            
                                        comp[i] = 1.0

                                if (mat[0] == "rectangle"):
                                    if ((mat[1] <= x_pos <= mat[2]) and (mat[3] <= y_pos <= mat[4])):
                                        for ii in range(0, i): # Erase previous composition
                                            comp[ii] = 0.0

                                        comp[i] = 1.0

                                if (mat[0] == "interface"):
                                    if (mat[1] == "below" and y_pos <= interface(x_pos)):
                                        for ii in range(0, i): # Erase previous composition
                                            comp[ii] = 0.0

                                        comp[i] = 1.0

                                    if (mat[1] == "above" and y_pos >= interface(x_pos)):
                                        for ii in range(0, i): # Erase previous composition
                                            comp[ii] = 0.0

                                        comp[i] = 1.0
                                i += 1

                        # Check whether the materials are assigned correctly
                        sum_comp = 0.0
                        if (len(comp) > 0):
                            for ii in range(0, len(comp)):
                                sum_comp += comp[ii]
                        
                            if (sum_comp != 1.0):
                                print("Sum of the compositions on a tracer does not equal 1!")
                                print("Exiting.")
                                MPI.comm_world.Abort()

                        # --- Create an ID that will be unique ---
                        unique_ID = tracer_no*size + rank

                        self.tracers.append([x_pos, y_pos, j, rank, 0.0, 0.0, 0.0, ocean_mat, surf_mat, y_pos, comp, [0.0, 0.0], 0.0, unique_ID])

                        # 0) x-position                     \
                        # 1) y-position                      \
                        # 2) cell where the tracer is         --- Essential for advection of tracers ---
                        # 3) rank where the tracer is        /

                        # --- Quantities to be advected ---
                        # 4) tau_xx
                        # 5) tau_xz
                        # 6) plastic strain
                        # 7) ocean material 
                        # 8) surface material
                        # 9) original depth of the tracer
                        # 10) composition [mat1, mat2, ...]
                        # 11) [melt fraction, delta_T]
                        # 12) tracer original (0) or added (1)
                        # 13) unique tracer ID

                        self.tracers_in_cells[j].append(tracer_no)
                        tracer_no +=1   

    def load_tracers(self):
        # infile = open("data_"+reload_name+"/tracers/step_"+str(reload_step)+".dat", "r") 
        infile = open("step_3000_new.dat", "r") 
        lines = infile.readlines() 
        tracer_no=0

        if (rank == 0):
            print("Reading tracers...")

        for j in range(self.mesh.num_cells()):
            self.tracers_in_cells.append([])

        header = True
        for line in lines:
            sline = line.split("\t")

            if (header==True):
                header = False
                continue

            x_pos = float(sline[0])
            y_pos = float(sline[1])

            try:
                # Test to see whether the point belongs to this rank
                test = self.Temp(Point(x_pos,y_pos))

                tau_xx = float(sline[2])
                tau_xz = float(sline[3])
                eps_plast = float(sline[4])
                ocean_mat = float(sline[5])
                surf_mat = float(sline[6])
                y_pos_orig = float(sline[7])
                comp = float(sline[8])
                mf = float(sline[9])
                origin = float(sline[10])

                for j in range(self.mesh.num_cells()):
                    if (Cell(self.mesh, j).contains(Point(x_pos, y_pos))):
                        self.tracers.append([x_pos, y_pos, j, rank, rank, tau_xx, tau_xz, eps_plast, ocean_mat, surf_mat, y_pos_orig, comp, mf, origin])
                        self.tracers_in_cells[j].append(tracer_no)
                        tracer_no += 1 
                        break # (once found, quit searching)
            except:
                pass

        MPI.barrier(comm)
        infile.close()

    def advect_tracers(self, v, v_mesh, dt):
        
        self.start_advection = time.time()
        tracers_process = len(self.tracers)
        vacancies_process = len(self.vacancy)

        tracers_total = MPI.sum(comm, tracers_process)
        vacancies_total = MPI.sum(comm, vacancies_process)

        if (rank == 0):
            print(" ")
            print("\tAdvecting", int(tracers_total - vacancies_total), "tracers...")

        self.vacancy.sort(reverse = False)
        self.moving_tracers = []

        MPI.barrier(comm)

        for j in range (0,len(self.tracers)):  
                if (len(self.tracers[j]) == 0): # Exclude vacant positions
                    continue

                px = self.tracers[j][0]
                py = self.tracers[j][1]
                v1 = v(Point(px,py))

                try:    
                    if (integration_method == "RK2"):
                        # --- Faster than RK4 and reasonable for most applications ---
                        v2 = v(Point(px + v1[0]*float(dt), py + v1[1]*float(dt)))

                        vx_total = (v1[0] + v2[0])/2.0
                        vy_total = (v1[1] + v2[1])/2.0

                    if (integration_method == "RK4"):
                        # --- Safe option ---
                        v2 = v(Point(px + v1[0]*float(dt)/2.0, py + v1[1]*float(dt)/2.0))
                        v3 = v(Point(px + v2[0]*float(dt)/2.0, py + v2[1]*float(dt)/2.0))
                        v4 = v(Point(px + v3[0]*float(dt),     py + v3[1]*float(dt)))

                        vx_total = (v1[0] + 2.0*v2[0] + 2.0*v3[0] + v4[0])/6.0
                        vy_total = (v1[1] + 2.0*v2[1] + 2.0*v3[1] + v4[1])/6.0

                    if (self.moving_mesh == True):
                        # If the final position is in the domain, will it be even after mesh  moves?
                        # Use the mesh velocity in the new position
                        vm = v_mesh(Point(px + vx_total*float(dt), py + vy_total*float(dt)))
                        
                        # Now the mesh move needs to be checked before updating the position.
                        # If the position is updated, but wouldn't be in the rank after mesh moved,
                        # it would go to the "except part" and would be advected two times !!!

                        # Test for new mesh position (must presede the assignment of new position)
                        v_test = v(Point(px + vx_total*float(dt), py + vy_total*float(dt) + 1.0*vm[1]*float(dt)))
                        v_test = v(Point(px + vx_total*float(dt), py + vy_total*float(dt) - 1.0*vm[1]*float(dt)))

                    # Only then update the tracer position (must be after the checking the mesh movement)
                    self.tracers[j][0] = px + vx_total*float(dt)
                    self.tracers[j][1] = py + vy_total*float(dt)

                    # Assign the rank at the final position of the tracer, not earlier!
                    v_test = v(Point(self.tracers[j][0], self.tracers[j][1]))
                    self.tracers[j][3] = rank

                except:

                    self.moving_tracers.append(self.tracers[j])
                    self.tracers_in_cells[self.tracers[j][2]].remove(j)
                    self.vacancy.append(j)
                    self.tracers[j] = []

        if (tracers_total > 0):
            print("\t---> rank", rank, "is advecting", "{:.1f}%.".format((tracers_process - vacancies_process) / (tracers_total - vacancies_total)*100.0))
        else:
            if (rank == 0):
                print("\t---> There are no tracers to advect.")

    def find_tracers(self, v, dt):

        self.tracers_to_find = []
        all_moving_tracers = comm.allgather(self.moving_tracers)

        for k in range(0, size):
            for j in range (0, len(all_moving_tracers[k])):
                px = all_moving_tracers[k][j][0]
                py = all_moving_tracers[k][j][1]

                count1 = 0
                count2 = 0
                count3 = 0
                count4 = 0
                count5 = 0

                if (integration_method == "RK2"):
                    try:
                        v1 = v(Point(px,py))
                    except:
                        v1 = [0.0, 0.0]
                        count1 += 1

                    v1_x    = MPI.sum(comm, v1[0])
                    v1_y    = MPI.sum(comm, v1[1])
                    count1  = MPI.sum(comm, count1)

                    try:
                        v2 = v(Point(px + v1_x*float(dt), py + v1_y*float(dt)))
                    except:
                        v2 = [0.0, 0.0]
                        count2 += 1
                    
                    v2_x    = MPI.sum(comm, v2[0])
                    v2_y    = MPI.sum(comm, v2[1])   
                    count2  = MPI.sum(comm, count2)        

                    all_moving_tracers[k][j][0] = px + (v1_x + v2_x)*float(dt)/2.0
                    all_moving_tracers[k][j][1] = py + (v1_y + v2_y)*float(dt)/2.0

                if (integration_method == "RK4"):
                    try:
                        v1 = v(Point(px,py))
                    except:
                        v1 = [0.0, 0.0]
                        count1 += 1

                    v1_x    = MPI.sum(comm, v1[0])
                    v1_y    = MPI.sum(comm, v1[1])
                    count1  = MPI.sum(comm,count1)

                    try:
                        v2 = v(Point(px + v1_x*float(dt)/2.0, py + v1_y*float(dt)/2.0))
                    except:
                        v2 = [0.0, 0.0]
                        count2 += 1
                    
                    v2_x    = MPI.sum(comm, v2[0])
                    v2_y    = MPI.sum(comm, v2[1])   
                    count2  = MPI.sum(comm, count2)  

                    try:
                        v3 = v(Point(px + v2_x*float(dt)/2.0, py + v2_y*float(dt)/2.0))
                    except:
                        v3 = [0.0, 0.0]
                        count3 += 1

                    v3_x    = MPI.sum(comm, v3[0])
                    v3_y    = MPI.sum(comm, v3[1])
                    count3  = MPI.sum(comm, count3)

                    try:
                        v4 = v(Point(px + v3_x*float(dt), py + v3_y*float(dt)))
                    except:
                        v4 = [0.0, 0.0]
                        count4 += 1

                    v4_x    = MPI.sum(comm, v4[0])             
                    v4_y    = MPI.sum(comm, v4[1])  
                    count4  = MPI.sum(comm, count4)       

                    all_moving_tracers[k][j][0] = px + (v1_x + 2.0*v2_x + 2.0*v3_x + v4_x)*float(dt)/6.0
                    all_moving_tracers[k][j][1] = py + (v1_y + 2.0*v2_y + 2.0*v3_y + v4_y)*float(dt)/6.0

                # Find out the rank at the final (!!) position of the tracer, not earlier!
                try:
                    v_test = v(Point(all_moving_tracers[k][j][0], all_moving_tracers[k][j][1]))
                    all_moving_tracers[k][j][3] = rank
                except:
                    all_moving_tracers[k][j][3] = 0
                    count5 += 1
                
                all_moving_tracers[k][j][3] = MPI.sum(comm, all_moving_tracers[k][j][3])
                count5 = MPI.sum(comm,count5) 

                not_found = False
                if (count1 == size or count2 == size or count3 == size or count4 == size or count5 == size):
                    # if count = size, the exeption has been performed by all processes, meaning that the tracer
                    # left the domain, whatever the reason
                    not_found = True

                if (rank == int(all_moving_tracers[k][j][3]) and not_found == False):
                    if (len(self.vacancy) > 0): # If possible, put them in the vacancies
                        self.tracers[self.vacancy[0]] = all_moving_tracers[k][j]
                        self.tracers_to_find.append(self.vacancy[0]) 
                        del self.vacancy[0]

                    else: # If no vacancy is available, put them at the end of the list
                        self.tracers.append(all_moving_tracers[k][j])
                        self.tracers_to_find.append(len(self.tracers) - 1) # The last added tracer
        
        self.end_advection = time.time()

        if (rank == 0):
            print("\tDone",'({:.2f} s).\n'.format(self.end_advection - self.start_advection))
            print(" ")

    def check_tracers(self):
        # --- Check number of tracers in cells ---
        TiC_sum = 0

        for j in range(self.mesh.num_cells()):
            TiC_sum += len(self.tracers_in_cells[j])

        TiC_sum     = MPI.sum(comm, TiC_sum)
        tracers_sum = MPI.sum(comm, len(self.tracers))
        vacancy_sum = MPI.sum(comm, len(self.vacancy))

        if (rank == 0):
            print("\tTracer lists length  =", int(tracers_sum))
            print("\tVacancy lists length =", int(vacancy_sum))
            print(" ")
            print("\tTotal number of traces     =", int(tracers_sum - vacancy_sum))
            print("\tNumber of tracers in cells =", int(TiC_sum))
            
    def delete_and_find(self):
        if (rank == 0):
            print("")
            print("\tSorting moving tracers...")
        
        # self.check_tracers()
        MPI.barrier(comm)

        # --- Searching for the tracers that moved to a neighboring cell ---
        for j in range(self.mesh.num_cells()):
            to_remove = []
            
            for i in range(0, len(self.tracers_in_cells[j])):
                tracer_no = self.tracers_in_cells[j][i]
                px = self.tracers[tracer_no][0]
                py = self.tracers[tracer_no][1]

                new_cell = self.mesh.bounding_box_tree().compute_first_entity_collision(Point(px, py))
                if (j == new_cell):
                    continue
                else:
                    to_remove.append(tracer_no)

                    if (new_cell <= self.mesh.num_cells()):
                        self.tracers_in_cells[new_cell].append(tracer_no)
                        self.tracers[tracer_no][2] = new_cell

                    else:
                        print("From cell to cell: tracer not found:", px, py)
                        self.vacancy.append(tracer_no)
                        self.tracers[tracer_no] = []

            # Tracers to remove
            for l in range(len(to_remove)):
                self.tracers_in_cells[j].remove(to_remove[l])
        
        # self.check_tracers()

        # --- Searching for the new cells of the tracers that left the process ---
        for i in range(0, len(self.tracers_to_find)):
            tracer_no = self.tracers_to_find[i]
            px = self.tracers[tracer_no][0]
            py = self.tracers[tracer_no][1]

            j = self.mesh.bounding_box_tree().compute_first_entity_collision(Point(px, py))
            if (j <= self.mesh.num_cells()):
                self.tracers_in_cells[j].append(tracer_no)
                self.tracers[tracer_no][2] = j
                
            else:
                print("From rank to rank: tracer not found:", px, py)
                
                self.vacancy.append(tracer_no)
                self.tracers[tracer_no] = []

        self.check_tracers()
        if (rank == 0):
            print("\tDone.")
            print("")

    def add_tracers(self, t, step):

        # --- Search the mesh for empty cells ---
        case = 0
        for j in range(self.mesh.num_cells()):
            if (len(self.tracers_in_cells[j]) == 0):
                case +=1

        case = MPI.sum(comm, case)

        if (rank == 0):
            file = open("data_" + self.name + "/empty_cells.dat", "a")
            file.write(("%.5E\t" + 2*"%d\t" + "\n")%(float(t)/time_units, step, case))
            file.close()
            
        for j in range(self.mesh.num_cells()):                
            centroid    = Cell(self.mesh, j).midpoint()
            cell_j      = Cell(self.mesh,j).get_vertex_coordinates()
            # --- Do not add tracers to the cells where the composition is mixed
            # or if the cell is allowed not to have tracers ---
            skip_cell = False

            # Composition case
            comp_cell = []
            for mat in self.composition:
                comp = mat(Point(centroid.x(), centroid.y()))
                comp_cell.append(comp)
                if (0.0 < comp < 1.0):
                    skip_cell = True

            # Empty cell case
            if (empty_cells_allowed == True):
                for empty in empty_cells_region:
                    if (empty[0] == "rectangle"):
                            if ((empty[1] <= centroid.x() <= empty[2]) and (empty[3] <= centroid.y() <= empty[4])):
                                skip_cell = True
            
            if (skip_cell == True):
                continue

            # --- If the cell is empty (and is not supposed to be), save the info into a file ---
            if (len(self.tracers_in_cells[j]) == 0):
                file = open("data_" + self.name + "/empty_cells.dat", "a")
                file.write((5*"\t" + 2*"%d\t" + 2*"%.3E\t" + "\n")%(j, rank, centroid.x(), centroid.y()))
                file.close()

            x_min   = min(cell_j[0], min(cell_j[2], cell_j[4]))
            x_max   = max(cell_j[0], max(cell_j[2], cell_j[4]))
            y_min   = min(cell_j[1], min(cell_j[3], cell_j[5]))
            y_max   = max(cell_j[1], max(cell_j[3], cell_j[5]))
            eps     = 1e-6

            stress_dev_xx_cell      = self.stress_dev_tensor(Point(centroid.x(),centroid.y()))[0]
            stress_dev_xz_cell      = self.stress_dev_tensor(Point(centroid.x(),centroid.y()))[1]
            eps_plast_cell          = self.plastic_strain(Point(centroid.x(),centroid.y()))
            om_cell                 = self.ocean_material(Point(centroid.x(),centroid.y()))
            sm_cell                 = self.surface_material(Point(centroid.x(),centroid.y()))
            xm_cell                 = self.xm(Point(centroid.x(),centroid.y()))
            topo_height             = self.h_top(Point(centroid.x(),centroid.y()))

            if (x_min >= eps and x_max <= length - eps): # Cells that are NOT on the vertical boundaries
                r_in    = Cell(self.mesh, j).inradius()
                deg2rad = np.pi/180.0
                lim     = int(sqrt(tracers_per_cell*4.0))

                # --- Number of tracers in the inner circle...
                n_1 = 0
                for i in range(0, len(self.tracers_in_cells[j])):
                    tracer_no = self.tracers_in_cells[j][i]
                    dx = abs(self.tracers[tracer_no][0] - centroid.x())
                    dy = abs(self.tracers[tracer_no][1] - centroid.y())
                    if (sqrt(dx**2 + dy**2) < 0.7*r_in):
                        n_1 += 1

                # ...and those outside the circle ---
                n_2 = len(self.tracers_in_cells[j]) - n_1

                # Add tracers only if 
                # 1) the number of tracers within the circle is below limit as well as the overall number
                # 2) the number outside the circle is lower than half of the limit 
                #    (protection from tracers that focus too much)

                if ((n_1 < 5 and len(self.tracers_in_cells[j]) < lim**2/4.0) or n_2 < lim**2/8.0):
                    for ii in range (0, 5):
                        tracer_angle = 90.0*(ii - 1)

                        if (ii==0):
                            # First tracer to the cell center
                            x_pos = centroid.x()
                            y_pos = centroid.y()
                        else:
                            # Other tracers along a circle
                            x_pos = centroid.x() + 0.5*r_in*sin(tracer_angle*deg2rad)
                            y_pos = centroid.y() + 0.5*r_in*cos(tracer_angle*deg2rad)

                            # --- Ocean material ---
                            if (y_pos < om_thickness or om_cell == 1.0):
                                om = 1.0
                            else:
                                om = 0.0

                            # --- Surface material ---
                            if (y_pos > ((height + topo_height) - sm_thickness) or sm_cell == 1.0):
                                sm = 1.0
                            else:
                                sm = 0.0

                            # --- Surface material and stripes ---
                            # The tracers that most likely falsely enter through the top boundary
                            if (y_pos > (height + topo_height) - sm_thickness):
                                original_y = height*self.height_fraction(Point(x_pos, y_pos))
                                origin = 2.0 # Added and distinguish from others

                            # --- Other tracers added in the rest of the domain, exclude them in both respects ---
                            else:
                                original_y = 0.0
                                origin = 1.0 # Treat it like added

                            self.tracers.append([x_pos, y_pos, j, rank, stress_dev_xx_cell, stress_dev_xz_cell,\
                                eps_plast_cell, om, sm, original_y, comp_cell, [xm_cell, 0.0], origin, 0.0])
                            
                            self.tracers_in_cells[j].append(len(self.tracers) - 1)  

            else: # Cells that are on the vertical boundaries
                adding = False

                # Among the left-boundary cells, find the leftmost tracer x-coordinate
                if (x_min < eps): 
                   
                    x_leftmost = length
                    for i in range(0, len(self.tracers_in_cells[j])):
                        tracer_no = self.tracers_in_cells[j][i]

                        if (self.tracers[tracer_no][0] < x_leftmost):
                            x_leftmost = self.tracers[tracer_no][0]

                    if (x_leftmost > centroid.x()):
                        adding = True
                        case = "left"

                # Among the right-boundary cells, find the rightmost tracer x-coordinate
                if (x_max > length - eps): 
                    x_rightmost = 0.0
                    for i in range(0, len(self.tracers_in_cells[j])):
                        tracer_no = self.tracers_in_cells[j][i]
                        if (self.tracers[tracer_no][0] > x_rightmost):
                            x_rightmost = self.tracers[tracer_no][0]

                    if (x_rightmost < centroid.x()):
                        adding = True
                        case = "right"

                # --- Add only when the last tracer is half the cell away ---
                if (adding == True):
                    for ii in range(0,10):

                        if (case == "left"):
                            x_pos = 1.0

                        elif (case == "right"):
                            x_pos = length - 1.0
                        
                        y_pos = y_min + ii/10.0*(y_max - y_min)

                        if (Cell(self.mesh, j).contains(Point(x_pos, y_pos)) == True):

                            # --- Ocean material ---
                            if (y_pos < om_thickness or om_cell == 1.0):
                                om = 1.0
                            else:
                                om = 0.0

                            # --- Surface material ---
                            if (y_pos > ((height + topo_height) - sm_thickness) or sm_cell == 1.0):
                                sm = 1.0
                            else:
                                sm = 0.0

                            original_y = height*self.height_fraction(Point(x_pos, y_pos))
                            origin = 0.0

                            self.tracers.append([x_pos, y_pos, j, rank, stress_dev_xx_cell, stress_dev_xz_cell,\
                                eps_plast_cell, om, sm, original_y, comp_cell, [xm_cell, 0.0], origin, 0.0])
                            
                            self.tracers_in_cells[j].append(len(self.tracers) - 1)