from dolfin import *
import numpy
from m_parameters import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

class MeshModule:
    def __init__(self):
        self.name = name

        if (reload == True):
            self.load_mesh()

        else:
            self.create_mesh()            
            self.refine_mesh()
            self.define_boundaries()

        self.x = SpatialCoordinate(self.mesh)
        self.ds = Measure("ds",subdomain_data = self.boundary_parts)

        self.e_z = Expression(("0.0", "1.0"), domain = self.mesh, degree = 1)
        self.e_x = Expression(("1.0", "0.0"), domain = self.mesh, degree = 1)
        self.normal = FacetNormal(self.mesh)

        self.bound_mesh = BoundaryMesh(self.mesh, 'exterior')

        if (BC_Stokes_problem[0][0] == "free_surface" or BC_Stokes_problem[1][0] == "free_surface"):
            self.moving_mesh = True
        else:
            self.moving_mesh = False

    def move_mesh(self, u_mesh):
        ALE.move(self.mesh, u_mesh)
        
        self.bound_mesh = BoundaryMesh(self.mesh, 'exterior') 

        self.mesh.init()
        self.mesh.bounding_box_tree().build(self.mesh)  
        self.mesh.init()

    def create_mesh(self):
        self.mesh = RectangleMesh(Point(0.0,0.0), Point(length, height), x_div, z_div, triangle_types)
    
    def refine_mesh(self):
        for i in range(0, len(refinement), 4):
            cell_markers = MeshFunction("bool", self.mesh, self.mesh.topology().dim())
            cell_markers.set_all(False)

            for cell in cells(self.mesh):
                mid_x = cell.midpoint().x()
                mid_y = cell.midpoint().y()
                
                if (mid_x > refinement[i] and mid_x < refinement[i+1] and mid_y > refinement[i+2] and mid_y < refinement[i+3]):
                    cell_markers[cell] = True
        
            self.mesh = refine(self.mesh, cell_markers)

    def define_boundaries(self):
        self.boundary_parts = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1, 0)

        top     = AutoSubDomain(lambda x: near(x[1], height))
        bottom  = AutoSubDomain(lambda x: near(x[1], 0.0))
        left    = AutoSubDomain(lambda x: near(x[0], 0.0))
        right   = AutoSubDomain(lambda x: near(x[0], length))

        top.mark(self.boundary_parts, 1)
        bottom.mark(self.boundary_parts, 2)
        left.mark(self.boundary_parts, 3)
        right.mark(self.boundary_parts, 4)

    def load_mesh(self):
        
        file_time = open(reload_name + "/HDF5/data_timestamp.dat")
        lines = file_time.readlines()

        if (reload_step == "last"):
            step = len(lines) - 2
        else:
            step = reload_step

        file_time.close()

        self.mesh = Mesh()
        self.mesh_file = HDF5File(comm, reload_name + "/HDF5/meshes/mesh_" + str(step) + ".h5", "r")
        self.mesh_file.read(self.mesh, "/mesh", True)
        
        self.boundary_parts = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1, 0)
        self.mesh_file.read(self.boundary_parts, "/subdomains")
        self.mesh_file.close()

def count_nodes(bound_mesh):
    i_top = 0
    i_bot = 0

    for vertex in vertices(bound_mesh):
        if (vertex.point()[1] == height):
            i_top += 1

        if (vertex.point()[1] == 0.0):
            i_bot += 1

    return i_top, i_bot

def detect_oscillations(mesh, bound_mesh, diff_coef, x_div_top):
    # --- Top boundary has to have equidistant nodes ---
    dx = length/x_div_top
    targets = []

    # --- 1/ Left to right, 2/ Right to left ---
    for jj in range(0, 2):

        # ---- Searching from left to right ---
        if (jj == 0):
            i_min, i_max = 1, x_div_top - 4
        
        # ---- Searching from right to left ---
        if (jj == 1):
            i_min, i_max = x_div_top, 4

        for i in range(i_min, i_max):
            h1 = 0
            h2 = 0
            h3 = 0
            h4 = 0
            h5 = 0

            for vertex in vertices(bound_mesh):
                if (vertex.point()[1] > height/2.0):
                    if (jj == 0):
                        if (vertex.point()[0] == i*dx):
                            h1 = vertex.point()[1]

                        if (vertex.point()[0] == (i + 1)*dx):
                            h2 = vertex.point()[1]

                        if (vertex.point()[0] == (i + 2)*dx):
                            h3 = vertex.point()[1]

                        if (vertex.point()[0] == (i + 3)*dx):
                            h4 = vertex.point()[1]

                        if (vertex.point()[0] == (i + 4)*dx):
                            h5 = vertex.point()[1]
                    else:
                        if (vertex.point()[0] == i*dx):
                            h1 = vertex.point()[1]

                        if (vertex.point()[0] == (i - 1)*dx):
                            h2 = vertex.point()[1]

                        if (vertex.point()[0] == (i - 2)*dx):
                            h3 = vertex.point()[1]

                        if (vertex.point()[0] == (i - 3)*dx):
                            h4 = vertex.point()[1]

                        if (vertex.point()[0] == (i - 4)*dx):
                            h5 = vertex.point()[1]

            h1 = MPI.max(comm, h1)
            h2 = MPI.max(comm, h2)
            h3 = MPI.max(comm, h3)
            h4 = MPI.max(comm, h4)
            h5 = MPI.max(comm, h5)

            MPI.barrier(comm)

            # --- In order to capture oscillations even in tilted areas ---
            global_angle = atan((h5-h1)/(4.0*dx))

            angle1 = atan((h2 - h1)/dx) - global_angle
            angle2 = atan((h3 - h2)/dx) - global_angle
            angle3 = atan((h4 - h3)/dx) - global_angle
            angle4 = atan((h5 - h4)/dx) - global_angle

            sign1 = numpy.sign(angle1)
            sign2 = numpy.sign(angle2)
            sign3 = numpy.sign(angle3)
            sign4 = numpy.sign(angle4)
            
            MPI.barrier(comm)

            h12 = sqrt((h2-h1)**2 + dx**2)*sin(angle1)
            h23 = sqrt((h3-h2)**2 + dx**2)*sin(angle2)
            h34 = sqrt((h4-h3)**2 + dx**2)*sin(angle3)
            h45 = sqrt((h5-h4)**2 + dx**2)*sin(angle4)

            # --- A threshold for the oscillations ---
            diff_min = max(abs(h12), max(abs(h23), max(abs(h34), abs(h45))))

            # --- Mark if there are oscillations or the element is too deformed ---
            if ((diff_min > 10.0 and sign1 == sign3 and sign2 == sign4 and sign1 == -sign2) or abs(h2 - h3) > dx):
                targets.append([(i+2)*dx, diff_min]) # approx in the middle of the oscillation
            
    # --- Remove duplicates from the targets ---
    targets = list(set(targets))

    ranks = []
    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()

        smooth = 0
        max_coef = 0
        for i in range(0, len(targets)):
            smooth += (topo_diff*targets[i][1]/10.0)*exp(-(centroid.x() - targets[i][0])**2/5e6)

            # --- Scales the diffusion coefficient according to the biggest oscillation ---
            if ((topo_diff*targets[i][1]/10.0) > max_coef):
                max_coef = (topo_diff*targets[i][1]/10.0)

        smooth = min(max_coef, smooth)
        ranks.append(smooth)
        
    ranks = numpy.array(ranks)                                                            
    diff_coef.vector().set_local(ranks)

def save_topography(bound_mesh, name, q_ice, q_water, step):
    if (BC_Stokes_problem[0][0] == "free_surface"):
        
        topography_top = []
        for vertex in vertices(bound_mesh):
            if (vertex.point()[1] > height*0.75):
                topography_top.append([vertex.point()[0], vertex.point()[1]])

        # --- Top topography ---
        all_topography_top = []
        gather_topography = comm.allgather(topography_top)
        for k in range(0,size):
            for j in range (0,len(gather_topography[k])):
                all_topography_top.append([gather_topography[k][j][0], gather_topography[k][j][1]])

        maximum_left = 0
        maximum_right = 0
        for j in range(len(all_topography_top)):
            if (all_topography_top[j][0] == 0.0 and all_topography_top[j][1] > maximum_left):
                maximum_left = all_topography_top[j][1]

            if (all_topography_top[j][0] == length and all_topography_top[j][1] > maximum_right):
                maximum_right = all_topography_top[j][1]

        all_topography_top.sort(key = lambda x: x[0])
        if (rank == 0):
            file = open("data_" + name + "/topography/step_" + str(step - 1) + "_top.dat","w")

            file.write((2*"%.7E\t" + "\n")%(0.0, maximum_left))
            for j in range (0, len(all_topography_top)):    
                if (0 < all_topography_top[j][0] < length):
                    file.write((2*"%.7E\t" + "\n")%(all_topography_top[j][0], all_topography_top[j][1]))
            file.write((2*"%.7E\t" + "\n")%(length, maximum_right))

            file.close()  

    if (BC_Stokes_problem[1][0] == "free_surface"):
        # --- Bottom topography ---
        topography_bot = []
        topography_bot_header = ["#x (m)\t\t", "y (m)\t\t"]
        for vertex in vertices(bound_mesh):
            if (vertex.point()[1] < height*0.25):
                x_pos =vertex.point()[0] 
                y_pos = vertex.point()[1]
                topography_bot.append([x_pos, y_pos])

                # --- Uncomment for the heat flux along the bottom bounrady ---
                # q_ice_vertex = q_ice(Point(x_pos, y_pos))[1]
                # q_water_vertex = q_water(Point(x_pos, y_pos))[1]
                # topography_bot[-1].append([q_ice_vertex, q_water_vertex])
                # topography_bot_header.append(["q_i (W/m2)\t", "q_w (W/m2)"])

        all_topography_bot = []
        gather_topography = comm.allgather(topography_bot)
        for k in range(0,size):
            for j in range (0,len(gather_topography[k])):
                all_topography_bot.append([gather_topography[k][j][0], gather_topography[k][j][1], gather_topography[k][j][2]])

        minimum_left = height
        minimum_right = height
        for j in range(len(all_topography_bot)):
            if (all_topography_bot[j][0] == 0.0 and all_topography_bot[j][1] < minimum_left):
                minimum_left = all_topography_bot[j][1]

            if (all_topography_bot[j][0] == length and all_topography_bot[j][1] < minimum_right):
                minimum_right = all_topography_bot[j][1]

        all_topography_bot.sort(key = lambda x: x[0])
        if (rank == 0):
            file = open("data_" + name + "/topography/step_" + str(step - 1) + "_bot.dat","w")
            file.write(((len(topography_bot_header[2]) + 2)*"%s\t" + "\n")%(topography_bot_header[0],
                                                                       topography_bot_header[1],
                                                                       *topography_bot_header[2]))

            # --- Left bottom point ----
            for j in range (len(all_topography_bot) - 1, -1, -1):  
                if (all_topography_bot[j][0] == 0 and all_topography_bot[j][1] == minimum_left):
                    write_topography_file(file, all_topography_bot, j)
                    break

            # --- Intermediate bottom points ---
            for j in range (0, len(all_topography_bot)):    
                if (0 < all_topography_bot[j][0] < length):
                    write_topography_file(file, all_topography_bot, j)
            
            # --- Right bottom point ---
            for j in range (0, len(all_topography_bot)):  
                if (all_topography_bot[j][0] == length and all_topography_bot[j][1] == minimum_right):
                    write_topography_file(file, all_topography_bot, j)
                    break
            
            file.close()

def write_topography_file(file, all_topography_bot, j):
    if (len(all_topography_bot[j][2]) == 0):
        file.write((2*"%.7E\t" + "\n")%(all_topography_bot[j][0], all_topography_bot[j][1]))
    else:
        file.write(((len(all_topography_bot[j][2]) + 2)*"%.7E\t" + "\n")%(all_topography_bot[j][0],
                                                                          all_topography_bot[j][1],
                                                                          *all_topography_bot[j][2]))