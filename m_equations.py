from dolfin import *
from m_parameters import *
from m_elements import *
from m_material_properties import *
from m_rheology import *
from m_boudnary_conditions import *
from m_timestep import *

class Equations:
    def __init__(self, MeshClass, ElemClass):

        self.sCG2 = ElemClass.sCG2
        self.sCG1 = ElemClass.sCG1
        self.vCG1 = ElemClass.vCG1
        self.sDG0 = ElemClass.sDG0
        self.tDG0 = ElemClass.tDG0
        self.V = ElemClass.V
        
        self.mesh = MeshClass.mesh
        self.boundary_parts = MeshClass.boundary_parts

        self.e_z = MeshClass.e_z
        self.normal = FacetNormal(MeshClass.mesh)


        # --- Trial and test functions ---
        self._Temp  = TrialFunction(self.sCG2)
        self.Temp_  = TestFunction(self.sCG2)
        self.Temp   = Function(self.sCG2)
        self.Temp_k = Function(self.sCG2)

        self.xm = Function(self.sDG0)
        self.visc = Function(self.sDG0)

        self.comp_0 = Function(self.sDG0)
        self.comp_1 = Function(self.sDG0)
        self.comp_2 = Function(self.sDG0)
        self.composition = [self.comp_0, self.comp_1, self.comp_2]

        self.number_of_tracers  = Function(self.sDG0)
        self.stress_dev_tensor  = Function(self.tDG0)
        self.plastic_strain     = Function(self.sDG0)
        self.ocean_material     = Function(self.sDG0)
        self.surface_material   = Function(self.sDG0)
        self.h_top              = Function(self.sDG0)
        self.mesh_ranks          = Function(self.sDG0)

        self.dt = Constant(1.0)

        self.top_length = Constant(1.0)
        self.unit_scalar = Function(self.sCG1)
        self.unit_scalar.assign(project(1.0, self.sCG1))
        self.height_fraction = Function(self.sCG1)
        self.height_fraction.assign(project(Expression("x[1]/height", height = height, degree=1), self.sCG1))

        self.q_top = Constant(1.0)
        self.q_cond_top = Constant(1.0)

        self._lambda = Constant(1.0)

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

        # --- Mesh displacement ---
        # Mesh displecement
        self._h2    = TrialFunction(self.sCG1)
        self.h2_    = TestFunction(self.sCG1)
        self.h2     = Function(self.sCG1)

        self.h2_top = Function(self.sCG1)     # Top Dirichlet BC for harmonic mesh displacement
        self.h2_bot = Function(self.sCG1)     # Bottom Dirichlet BC for harmonic mesh displacement


        

        if (stokes_elements == "Mini"):
            _pv = TrialFunction(self.V)
            self._p, _vl, _vb = split(_pv)
            self._v = _vl + _vb 

            pv_ = TestFunction(self.V)
            self.p_, vl_, vb_ = split(pv_)
            self.v_ = vl_ + vb_ 

            self.pv = Function(self.V)
            self.p, vl, vb = split(self.pv)
            self.v = vl + vb 

            self.v_k = Function(self.vCG1)
            self.v_mesh = Function(self.vCG1)
            self.u_mesh = Function(self.vCG1)

        if (stokes_elements == "TH"):
            _pv = TrialFunction(self.V)
            _p, _v = split(_pv)

            pv_ = TestFunction(self.V)
            p_, v_ = split(pv_)

            pv = Function(self.V)
            p, v = split(pv)

        top_left = Point_Fixed_Pressure()
        self.bc_temp = apply_temperature_BC(self.sCG2, self.boundary_parts)
        self.bc_stokes = apply_velocity_BC(self.V, self.boundary_parts, top_left)
       
        self.stationary_heat_eq()
        self.PerturbProfile()
        self.heat_eq()

        self.Stokes_eq()

        self.Equation_top_free_surface(self)
        self.Equation_bottom_free_surface(self)
        self.Equation_displacement_distribution(self)

    
    def stationary_heat_eq(self):
        # --- Boundary conditions ---
        bc_T_top = DirichletBC(self.sCG2, T_top, self.boundary_parts,1)
        bc_T_bot = DirichletBC(self.sCG2, T_bot, self.boundary_parts,2)
        self.bc_temp = [bc_T_bot, bc_T_top]

        # --- Equation ---
        eq_energy_ini      = dot(nabla_grad(self._Temp),nabla_grad(self.Temp_))*dx
        problem_energy_ini = LinearVariationalProblem(lhs(eq_energy_ini),rhs(eq_energy_ini), self.Temp, self.bc_temp)
        self.solver_energy_ini  = LinearVariationalSolver(problem_energy_ini)

    def PerturbProfile(self):
        self.Temp_k.assign(project(Expression("Temp - pa*sin(pi*x[1]/h)*cos(2*pf*pi*x[0]/l)", Temp = self.Temp, l=length, h=height, pa=perturb_ampl, pf=perturb_freq, degree=1), self.sCG2)) 
        self.Temp.assign(self.Temp_k) 

    def heat_eq(self):
        eq_energy = (rho_s*cp(self.Temp_k, self.composition)*(self._Temp - self.Temp_k)/self.dt*self.Temp_
                     
                + 0.5*(rho_s*cp(self.Temp_k, self.composition)*dot(self.v_k - self.v_mesh,nabla_grad(self.Temp))*self.Temp_ 
                       + k(self.Temp_k, self.composition)*dot(nabla_grad(self._Temp), nabla_grad(self.Temp_)))

                + 0.5*(rho_s*cp(self.Temp_k, self.composition)*dot(self.v_k - self.v_mesh,nabla_grad(self.Temp_k))*self.Temp_ 
                       + k(self.Temp_k, self.composition)*dot(nabla_grad(self.Temp_k),nabla_grad(self.Temp_))))*dx
        
        problem_energy = LinearVariationalProblem(lhs(eq_energy),rhs(eq_energy), self.Temp, self.bc_temp)
        self.solver_energy  = LinearVariationalSolver(problem_energy)

    def Stokes_eq(self):
        eq_cont = div(self._v)*self.p_*dx

        eq_momentum = (self._p*div(self.v_)
                - density(self.Temp, self.composition, self.xm)*g*dot(self.e_z, self.v_)
                - self.visc*inner(2.0*sym(nabla_grad(self._v)), nabla_grad(self.v_).T))*dx
        
        # --- Top surface "drunken sailor" stabilization ---
        if (BC_vel_top == "free_surface"): 
            eq_momentum += -density(self.Temp, self.xm, self.composition)*g*self._lambda*self.dt*dot(self._v, self.normal)*dot(self.ez, self.v_)*self.ds(1)

        # --- Hydrostatic pressure and bottom surface "drunken sailor"---
        if (BC_vel_bot == "free_surface" or phase_transition == True): 
            eq_momentum +=  dot(self.normal, self.v_)*(self.h_bot*rho_l - height*rho_s)*g*self.ds(2)\
                        + (rho_l - rho_s)*g*self._lambda*self.dt*dot(self._v, self.normal)*dot(self.ez, self.v_)*self.ds(2)

        problem_stokes = LinearVariationalProblem(lhs(eq_cont + eq_momentum),rhs(eq_cont + eq_momentum), self.pv, self.bc_stokes)
        self.solver_stokes = LinearVariationalSolver(problem_stokes)

    def Equation_top_free_surface(self):
            # --- Top free surface equation ---
        eq_fs_top = (self._h_top.dx(1)*self.h_top_.dx(1))*dx\
            - (self._h_top.dx(1)*self.h_top_)*self.ds(1)\
            - ((self._h_top - self.h_top_k + self.dt*(self._h_top.dx(0)*self.v_k[0] - self.v_k[1]))\
            * (self.h_top_.dx(1) - self.h_top_/(self.gamma*self.h_e_top)))*self.ds(1)

        bc_fs_top = []

        free_surface_top = LinearVariationalProblem(lhs(eq_fs_top), rhs(eq_fs_top), self.h_top, bc_fs_top)
        self.solver_free_surface_top = LinearVariationalSolver(free_surface_top)

    def Equation_bottom_free_surface(self):

        # --- Bottom free surface equation ---
        eq_fs_bot = (self._h_bot.dx(1))*(self.h_bot_.dx(1))*dx\
            - (self._h_bot.dx(1))*self.h_bot_*self.ds(2)\
            - (self._h_bot - self.h_bot_k + self.dt*(self._h_bot.dx(0)*(self.v_k[0] - self.u[0]) - (self.v_k[1] - self.u[1])))\
            * (self.h_bot_.dx(1) - self.h_bot_/(self.gamma*self.h_e_bot))*self.ds(2)
                                                                                    
        bc_fs_bot = []

        free_surface_bot = LinearVariationalProblem(lhs(eq_fs_bot), rhs(eq_fs_bot), self.h_bot, bc_fs_bot)
        self.solver_free_surface_bot = LinearVariationalSolver(free_surface_bot)

    # def Equation_topography_diffusion(self):

    #     # --- Topography diffusion equation ---
    #     eq_surf_diff = (_h_top-h_top1)/dt*h_top_*dx + diff_coef*dot(nabla_grad(_h_top),nabla_grad(h_top_))*dx
    #     bc_sd = []

    #     problem_surf_diff = LinearVariationalProblem(lhs(eq_surf_diff),rhs(eq_surf_diff),h_top,bc_sd)
    #     solver_surf_diff = LinearVariationalSolver(problem_surf_diff)

    def Equation_displacement_distribution(self):
        # --- Harmonics mesh displacement distribution ---
        eq_surf_move = dot(nabla_grad(self._h2),nabla_grad(self.h2_))*dx

        bc_sm = [DirichletBC(self.sCG1, self.h2_top, self.boundary_parts, 1),\
                 DirichletBC(self.sCG1, self.h2_bot, self.boundary_parts, 2)]

        problem_surf_move = LinearVariationalProblem(lhs(eq_surf_move), rhs(eq_surf_move), self.h2, bc_sm)
        self.solver_surf_move = LinearVariationalSolver(problem_surf_move)

    # ---------------------------------------------------------------
    # DEF call solvers 
    # ---------------------------------------------------------------
    def solve_heat_equation(self):
        self.solver_energy.solve()
        self.Temp_k.assign(self.Temp)

    def solve_Stokes_problem(self):
        self.visc.assign(project(1.0, self.sDG0))
        
        self.solver_stokes.solve()
        if (stokes_elements == "Mini"):
            self.p_out, self.v_out, self.vb = self.pv.split(True)
        if (stokes_elements == "TH"):
            self.p_out, self.v_out = self.pv.split(True)

        self.v_k.assign(self.v_out)

        self.dt.assign(time_step_domain(self.mesh, self.v_k, H_max, self.composition))

    def Solver_Topography_evolution(self):
        if (BC_vel_top == "free_surface"):
            self.solver_free_surface_top.solve()

        if (BC_vel_bot == "free_surface"):
            self.solver_free_surface_bot.solve()

        # --- Mesh displacement solver ---
        self.h2_top.assign(self.h_top - self.h_top_k)
        self.h2_bot.assign(self.h_bot - self.h_bot_k)

        self.h_top_k.assign(self.h_top)
        self.h_bot_k.assign(self.h_bot)

        self.solver_surf_move.solve()
        # self.u_mesh.assign(project(Expression(("0.0", "h2"), h2 = self.h2, degree = 1), self.vCG1))
        self.u_mesh.assign(project(self.h2*self.e_z, self.vCG1))
        self.v_mesh.assign(project(self.u_mesh/self.dt, self.vCG1))