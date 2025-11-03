from dolfin import *
from m_parameters import *

class Elements:
    
    def __init__(self, MeshClass):
        self.mesh = MeshClass.mesh
        
        sDG0_elem = FiniteElement("DG", self.mesh.ufl_cell(), 0) 
        self.sDG0 = FunctionSpace(self.mesh, sDG0_elem)

        tDG0_elem = TensorElement("DG", self.mesh.ufl_cell(), 0) 
        self.tDG0 = FunctionSpace(self.mesh, tDG0_elem)

        vCG1_elem = VectorElement("Lagrange", self.mesh.ufl_cell(), 1)  
        self.vCG1 = FunctionSpace(self.mesh, vCG1_elem)

        vCG2_elem = VectorElement("Lagrange", self.mesh.ufl_cell(), 2)
        self.vCG2 = FunctionSpace(self.mesh, vCG2_elem)

        sCG1_elem = FiniteElement("Lagrange", self.mesh.ufl_cell(), 1)
        self.sCG1 = FunctionSpace(self.mesh, sCG1_elem)

        sCG2_elem = FiniteElement("Lagrange", self.mesh.ufl_cell(), 2)
        self.sCG2 = FunctionSpace(self.mesh, sCG2_elem)

        if (stokes_elements == "Mini"):
            vB_elem  = VectorElement("B", self.mesh.ufl_cell(), self.mesh.topology().dim() + 1)
            self.V = FunctionSpace(self.mesh, MixedElement([sCG1_elem, vCG1_elem, vB_elem]))
            self.stokes_space = self.sCG1
            
        if (stokes_elements == "TH"):
            self.V = FunctionSpace(self.mesh, MixedElement([sCG1_elem, vCG2_elem]))
            self.stokes_space = self.sCG2

        self._Temp  = TrialFunction(self.sCG2)
        self.Temp_  = TestFunction(self.sCG2)
        self.Temp   = Function(self.sCG2)
        self.Temp_k = Function(self.sCG2)

        self.xm = Function(self.sDG0)
        self.melting_rate = Function(self.sDG0)
        self.xm_k = Function(self.sDG0)
        self.heating = Function(self.sDG0)
        self.heating_shear = Function(self.sDG0)
        self.visc = Function(self.sDG0)
        self.eta_v = Function(self.sDG0)
        self.log10_visc = Function(self.sDG0)
        self.delta_T = Function(self.sDG0)

        self.plastic_strain     = Function(self.sDG0)
        self.yield_stress     = Function(self.sDG0)
        self.yield_function     = Function(self.sDG0)
        

        self.comp_0 = Function(self.sDG0)
        self.comp_1 = Function(self.sDG0)
        self.comp_2 = Function(self.sDG0)
        self.composition = [self.comp_0, self.comp_1, self.comp_2]

        self.number_of_tracers  = Function(self.sDG0)
        self.stress_dev_tensor  = Function(self.tDG0)
        self.plastic_strain     = Function(self.sDG0)
        self.ocean_material     = Function(self.sDG0)
        self.surface_material   = Function(self.sDG0)
        self.mesh_ranks         = Function(self.sDG0)
        self.density            = Function(self.sDG0)
        self.grain_size         = Function(self.sDG0)

        self.stress_dev_inv      = Function(self.sDG0)
        self.stress_dev_inv_k     = Function(self.sDG0)
        self.strain_rate_inv     = Function(self.sDG0)
        
        self.iteration_error     = Function(self.sCG1)

        # --- Elasticity --
        self.stress_dev_xx     = Function(self.sDG0)
        self.stress_dev_xz      = Function(self.sDG0)
        self.z_function          = Function(self.sDG0)
        self.shear_modulus      = Function(self.sDG0)

        self.strain_rate_tensor = Function(self.tDG0)
        self.stress_dev_tensor = Function(self.tDG0)
        self.stress_dev_tensor_k = Function(self.tDG0)
        self.vorticity_tensor = Function(self.tDG0)

        self.dt = Constant(1.0)

        self.top_length = Constant(1.0)
        self.unit_scalar = Function(self.sCG1)
        self.unit_scalar.assign(project(1.0, self.sCG1))
        self.height_fraction = Function(self.sCG1)
        self.height_fraction.assign(project(Expression("x[1]/height", height = height, degree=1), self.sCG1))

        self.q_top = Constant(1.0)
        self.q_cond_top = Constant(1.0)

        self._lambda = Constant(1.0)

        self.v_k = Function(self.vCG1)
        self.u = Function(self.vCG1)
        self.p_k      = Function(self.sCG1)
        self.u_mesh = Function(self.vCG1)
        self.v_mesh = Function(self.vCG1)

        self.q_ice = Function(self.vCG1)
        self.q_water = Function(self.vCG1)
        self.n_bot = Function(self.vCG1)

        # --- Top and bottom boundary topography ---
        self._h_top     = TrialFunction(self.sCG1)
        self.h_top_     = TestFunction(self.sCG1)
        self.h_top      = Function(self.sCG1)
        
        self.h_top_k    = Function(self.sCG1)
        self.h_top1     = Function(self.sCG1)

        self._h_bot     = TrialFunction(self.sCG1)
        self.h_bot_     = TestFunction(self.sCG1)
        self.h_bot      = Function(self.sCG1)
        self.h_bot_k    = Function(self.sCG1)
        self.h_bot_pr    = Function(self.sCG1)

        # --- Mesh displacement ---
        # Mesh displecement
        self._h2    = TrialFunction(self.sCG1)
        self.h2_    = TestFunction(self.sCG1)
        self.h2     = Function(self.sCG1)

        self.h2_top = Function(self.sCG1)     # Top Dirichlet BC for harmonic mesh displacement
        self.h2_bot = Function(self.sCG1)     # Bottom Dirichlet BC for harmonic mesh displacement

        # --- iSALE ---
        self.Temp_iSALE     = Function(self.sDG0)
        self.h_top_iSALE    = Function(self.sDG0)
        self.density_iSALE  = Function(self.sDG0)