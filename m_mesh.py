from dolfin import *
import numpy
from m_parameters import *

comm = MPI.comm_world
rank = MPI.rank(comm)
size = MPI.size(comm)

class MeshModule:
    def __init__(self):
        if (loading_mesh == True):
            self.Load_Mesh()
        else:
            self.create_mesh()            
            self.refine_mesh()
            self.define_boundaries()

        self.save_mesh()

        self.x = SpatialCoordinate(self.mesh)
        self.ds = Measure("ds",subdomain_data = self.boundary_parts)

        self.e_z = Expression(("0.0", "1.0"), domain = self.mesh, degree = 1)
        self.e_x = Expression(("1.0", "0.0"), domain = self.mesh, degree = 1)
        self.normal = FacetNormal(self.mesh)

        self.bound_mesh = BoundaryMesh(self.mesh, 'exterior')

        if (BC_vel_top == "free_surface" or BC_vel_bot == "free_surface"):
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
    
    def save_mesh(self):
        mesh_file = XDMFFile(comm, "data_"+name+"/HDF5/mesh.xdmf")
        mesh_file.write(self.mesh)
        mesh_file.close()

        file = File("data_"+name+"/HDF5/subdomains.pvd")
        file.write(self.boundary_parts)

    def Save_Mesh(self, step_output):
        mesh_file = HDF5File(comm, "data_" + name + "/HDF5/meshes/mesh_" + str(int(step_output - 1)) + ".h5", "w")
        mesh_file.write(self.mesh, "/mesh")
        mesh_file.close()

    def Load_Mesh(self):
        self.mesh = Mesh()
        self.mesh_file = HDF5File(comm, "data_" + name + "/HDF5/meshes/mesh_" + str(int(reload_HDF5)) + ".h5", "r")
        self.mesh_file.read(mesh, "/mesh", True)
        self.mesh_file.close()

    def define_boundaries(self):
        self.boundary_parts = MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1, 0)

        top     = AutoSubDomain(lambda x: near(x[1], height))
        bottom  = AutoSubDomain(lambda x: near(x[1], 0.0))
        left    = AutoSubDomain(lambda x: near(x[0], 0.0))
        right   = AutoSubDomain(lambda x: near(x[0], length))

        top.mark(self.boundary_parts, 1)
        bottom.mark(self.boundary_parts, 2)
        left.mark(self.boundary_parts, 3)
        right.mark(self.boundary_parts, 4)

def output_timing(step, time_out):
    
    if (step == 1 or (output_frequency[0] == "steps" and step % output_frequency[1] == 0)\
        or (output_frequency[0] == "time" and float(time_out) >= output_frequency[1])):
        
        value = True
    else:
        value = False

    return value

def count_nodes(bound_mesh):
    i_top = 0
    i_bot = 0

    for vertex in vertices(bound_mesh):
        if (vertex.point()[1] == height):
            i_top += 1

        if (vertex.point()[1] == 0.0):
            i_top += 1

    return i_top, i_bot

def detect_oscillations(mesh, bound_mesh, diff_coef, x_div_top):
    # Top boundary has to have equidistant nodes
    dx = length/x_div_top
    nmax = x_div_top-4

    targets = []

    for i in range(1,nmax):

        h1 = 0
        h2 = 0
        h3 = 0
        h4 = 0
        h5 = 0

        for vertex in vertices(bound_mesh):
            if (vertex.point()[1]>height/2.0):
                if (vertex.point()[0] == i*dx):
                    h1 = vertex.point()[1]
                if (vertex.point()[0] == (i+1)*dx):
                    h2 = vertex.point()[1]
                if (vertex.point()[0] == (i+2)*dx):
                    h3 = vertex.point()[1]
                if (vertex.point()[0] == (i+3)*dx):
                    h4 = vertex.point()[1]
                if (vertex.point()[0] == (i+4)*dx):
                    h5 = vertex.point()[1]

        h1 = MPI.max(comm, h1)
        h2 = MPI.max(comm, h2)
        h3 = MPI.max(comm, h3)
        h4 = MPI.max(comm, h4)
        h5 = MPI.max(comm, h5)

        MPI.barrier(comm)

        # In order to capture oscillations even in tilted areas
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

        # A threshold for the oscillations
        diff_min = max(abs(h12),max(abs(h23),max(abs(h34),abs(h45))))

        # Mark if there are oscillations or the element is too deformed
        if ((diff_min>10.0 and sign1==sign3 and sign2==sign4 and sign1==-sign2) or abs(h2-h3)>dx):
            targets.append([(i+2)*dx,diff_min]) # approx in the middle of the oscillation
            
    ranks = []
    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()

        smooth = 0
        max_coef = 0
        for i in range(0,len(targets)):
            smooth += (topo_diff*targets[i][1]/10.0)*exp(-(centroid.x()-targets[i][0])**2/5e6)

            # Scales the diffusion coefficient according to the biggest oscillation
            if ((topo_diff*targets[i][1]/10.0) > max_coef):
                max_coef = (topo_diff*targets[i][1]/10.0)

        smooth = min(max_coef, smooth)
        ranks.append(smooth)
        
    ranks = numpy.array(ranks)                                                            
    diff_coef.vector().set_local(ranks)