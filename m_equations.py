from dolfin import *
from m_parameters import *
from m_elements import *
from m_material_properties import *
from m_rheology import *
from m_boundary_conditions import *
from m_timestep import *
from m_mesh import *
from m_interpolation import *

class Equations:
    def __init__(self, MeshClass, ElemClass, TracersClass, MeltingClass, FilesClass, t):

        self.sCG2   = ElemClass.sCG2
        self.sCG1   = ElemClass.sCG1
        self.vCG1   = ElemClass.vCG1
        self.sDG0   = ElemClass.sDG0
        self.tDG0   = ElemClass.tDG0
        self.V      = ElemClass.V

        self.stokes_space = ElemClass.stokes_space

        self.t = t

        if (TracersClass.use_tracers == True):
            self.tracers = TracersClass.tracers
            self.tracers_in_cells = TracersClass.tracers_in_cells
            self.only_melt_tracers = TracersClass.only_melt_tracers
            self.vacancy = TracersClass.vacancy

            self.add_melt_tracers       = MeltingClass.add_melt_tracers
            self.delete_melt_tracers    = MeltingClass.delete_melt_tracers
            self.compute_melting    = MeltingClass.compute_melting

        self.mesh   = MeshClass.mesh
        self.ds     = MeshClass.ds
        self.boundary_parts = MeshClass.boundary_parts

        if (viscosity_type == "GK_2001" or plasticity == True):
            self.nonlinear_rheology = True
        else:
            self.nonlinear_rheology = False

        self.viscosity_from = "strain_rate"
        if (elasticity == True):
            self.viscosity_from = "stress"

        self.e_z = MeshClass.e_z
        self.e_x = MeshClass.e_x
        self.normal = FacetNormal(MeshClass.mesh)

        # --- Trial and test functions ---
        self._Temp  = ElemClass._Temp
        self.Temp_  = ElemClass.Temp_
        self.Temp   = ElemClass.Temp
        self.delta_T   = ElemClass.delta_T
        self.Temp_k = ElemClass.Temp_k

        self.xm     = ElemClass.xm
        self.xm_k   = ElemClass.xm_k
        self.melting_rate   = ElemClass.melting_rate
        self.heating = ElemClass.heating
        self.visc       = ElemClass.visc
        self.eta_v       = ElemClass.eta_v
        self.log10_visc = ElemClass.log10_visc
        self.delta_T = ElemClass.delta_T

        self.stress_dev_inv      = ElemClass.stress_dev_inv
        self.stress_dev_inv_k      = ElemClass.stress_dev_inv_k
        self.strain_rate_inv     = ElemClass.strain_rate_inv

        self.stress_dev_xx      = ElemClass.stress_dev_xx
        self.stress_dev_xz      = ElemClass.stress_dev_xz
        self.z_function         = ElemClass.z_function
        self.shear_modulus         = ElemClass.shear_modulus

        self.strain_rate_tensor     = ElemClass.strain_rate_tensor
        self.stress_dev_tensor      = ElemClass.stress_dev_tensor
        self.vorticity_tensor       = ElemClass.vorticity_tensor

        self.u = ElemClass.u
        self.n_bot = ElemClass.n_bot
        self.comp_0 = ElemClass.comp_0
        self.comp_1 = ElemClass.comp_1
        self.comp_2 = ElemClass.comp_2
        self.composition = [self.comp_0, self.comp_1, self.comp_2]

        self.q_ice = ElemClass.q_ice
        self.q_water = ElemClass.q_water

        self.density     = ElemClass.density

        self.plastic_strain     = ElemClass.plastic_strain
        self.yield_stress       = ElemClass.yield_stress
        self.yield_function     = ElemClass.yield_function
        self.yield_function.assign(project(Expression("0.0", degree=1), self.sDG0))


        self.number_of_tracers  = Function(self.sDG0)
        
        self.ocean_material     = Function(self.sDG0)
        self.surface_material   = Function(self.sDG0)
        self.mesh_ranks          = Function(self.sDG0)

        if (time_step_strategy == "constant"):
            self.dt = Constant(dt_const)
        else:
            self.dt = Constant(0.1*kyr)

        self.h_bot_aver = Constant(0.0)
        self.q_ice_aver = Constant(0.0)
        self.q_ice_aver_ini = Constant(0.0)
        self.q_water_aver = Constant(0.0)

        self.sr_min = Constant(0.0)

        self.gamma = 0.005 
        self.h_e_top = length/x_div
        self.h_e_bot = length/x_div

        self.top_length = Constant(1.0)
        self.unit_scalar = Function(self.sCG1)
        self.unit_scalar.assign(project(1.0, self.sCG1))
        self.height_fraction = Function(self.sCG1)
        self.height_fraction.assign(project(Expression("x[1]/height", height = height, degree=1), self.sCG1))

        self.q_top = Constant(1.0)
        self.q_cond_top = Constant(1.0)

        self._lambda = Constant(1.0)

        self.v_aver = Constant(1.0)
        self.v_top_aver = Constant(1.0)

        # --- Top and bottom boundary topography ---        
        # self.h_top1     = Function(self.sCG1)

        self._h_top     = ElemClass._h_top
        self.h_top_     = ElemClass.h_top_
        self.h_top      = ElemClass.h_top
        self.h_top_k    = ElemClass.h_top_k

        self._h_bot     = ElemClass._h_bot
        self.h_bot_     = ElemClass.h_bot_
        self.h_bot      = ElemClass.h_bot
        self.h_bot_k    = ElemClass.h_bot_k

        self.iteration_error     = ElemClass.iteration_error

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

        if (stokes_elements == "TH"):
            _pv = TrialFunction(self.V)
            self._p, self._v = split(_pv)

            pv_ = TestFunction(self.V)
            self.p_, self.v_ = split(pv_)

            self.pv = Function(self.V)
            self.p, self.v = split(self.pv)

        self.name   = FilesClass.name
        self.v_k = ElemClass.v_k
        self.p_k = ElemClass.p_k
        self.v_mesh = ElemClass.v_mesh
        self.u_mesh = ElemClass.u_mesh

        top_left = Point_Fixed_Pressure()
        self.bc_temp    = apply_temperature_BC(self.sCG2, self.boundary_parts)
        self.bc_stokes  = apply_velocity_BC(self.V, self.boundary_parts, top_left)

        self.equation_heat_ini()
        self.thermal_perturbation()
        self.equation_heat()

        self.equation_Stokes()

        self.equation_free_surface_top()
        self.equation_free_surface_bottom()
        self.equation_mesh_displacement()
    
    def equation_heat_ini(self):
        """ 
        :var: Defines the weak form of the stationary heat transfer equation.

        .. math::
        
            0 = \\int_\\Omega k\\ (\\nabla T \\cdot \\nabla \\hat{T})\\ \\textrm{d}\\Omega
        
        If ``BC_heat_transfer[0][0] = "radiation"``, the corresponding boundary condition is added

        .. math::
        
            - \\int_{\\Gamma_{\\rm top}} I(1 - A)\\hat{T}\\ + T^4\\varepsilon\\sigma_{SB}\\hat{T}\\ \\textrm{d}\\Gamma_{\\rm top}

        where

        * :math:`I` is the insolation at the body's surface
        * :math:`A` is the albedo of the body's surface
        * :math:`\\varepsilon` is the body's surface emmisivity
        * :math:`\\sigma_{SB}` is the Stefan-Boltzmann constant

        .. note::

            If ``BC_heat_transfer[0][0] = "radiation"``, nonlinear solver will be automatically used.

        """
        # --- Equation ---
        # if (nonlinear_heat_equation == True):
        eq_energy_ini =  k(self.Temp, self.composition)*inner(nabla_grad(self.Temp), nabla_grad(self.Temp_))*dx

        if (BC_heat_transfer[0][0] == "radiation"):
            eq_energy_ini += - (1.0 - albedo)*insolation*self.Temp_*self.ds(1)\
                        + pow(self.Temp, 4)*emis*SB_const*self.Temp_*self.ds(1)
            
        if (initial_tidal_dissipation == True):
            eq_energy_ini += - tidal_heating(eta_eff(self.Temp, self.xm))*self.Temp_*dx
            
        J_eq_energy_ini = derivative(eq_energy_ini, self.Temp)
        problem_energy_ini = NonlinearVariationalProblem(eq_energy_ini, self.Temp, self.bc_temp, J_eq_energy_ini)
        self.solver_energy_ini  = NonlinearVariationalSolver(problem_energy_ini)

        prm = self.solver_energy_ini.parameters
        prm["newton_solver"]["absolute_tolerance"] = 1E-8
        prm["newton_solver"]["relative_tolerance"] = 1E-7
        prm["newton_solver"]["maximum_iterations"] = 25
        prm["newton_solver"]["relaxation_parameter"] = 1.0
        prm["newton_solver"]["linear_solver"] = "bicgstab"

        # else:
        #     eq_energy_ini  = dot(nabla_grad(self._Temp),nabla_grad(self.Temp_))*dx
        #     if (initial_tidal_dissipation == True):
        #         eq_energy_ini += -tidal_heating(self.visc)*self.Temp_*dx

        #     problem_energy_ini = LinearVariationalProblem(lhs(eq_energy_ini),rhs(eq_energy_ini), self.Temp, self.bc_temp)
        #     self.solver_energy_ini  = LinearVariationalSolver(problem_energy_ini)

    def thermal_perturbation(self):
        """ 
        :var: Perturbs the initial temperature profile.
        
        :param perturb_ampl: amplitude of the thermal perturbation (:math:`A`\ )
        :param perturb_freq: spatial frequency of the thermal perturbation (:math:`f`\ )

        :returns:
            .. math::
            
                T^\\prime = T - A\\sin\\left(\\frac{\\pi y}{h}\\right) \\cdot\\cos\\left(\\frac{2\\pi fx}{l}\\right)

        """
        self.Temp_k.assign(project(Expression("Temp - pa*sin(pi*x[1]/h)*cos(2*pf*pi*x[0]/l)", Temp = self.Temp, l=length, h=height, pa=perturb_ampl, pf=perturb_freq, degree=1), self.sCG2)) 
        self.Temp.assign(self.Temp_k) 

    def equation_heat(self):
        """ 
        :var: Defines the weak form of the time-dependent heat transfer equation.

        :returns:

        .. math::
        
            \\begin{align}
            0 = &\\int_\\Omega \\rho_s c_p \\frac{T - T^k}{\\Delta t}\\hat{T}\\\\
            &+ 0.5[\\rho_s c_p((\\boldsymbol{v} - \\boldsymbol{v}_{\\rm mesh}) \\cdot \\nabla T)\\ \\hat{T} + k\\ (\\nabla T \\cdot \\nabla \\hat{T})]\\\\
            &+ 0.5[\\rho_s c_p((\\boldsymbol{v} - \\boldsymbol{v}_{\\rm mesh}) \\cdot \\nabla T^k)\\ \\hat{T} + k\\ (\\nabla T^k \\cdot \\nabla \\hat{T})]\\ \\textrm{d}\\Omega\\\\
            \\end{align}
        
        where

        * :math:`T` is the temperature trial function
        * :math:`T^k` is the temperature function from the previous time step
        * :math:`\\hat{T}` is the temperature test function
        * :math:`\\boldsymbol{v}` is the velocity
        * :math:`\\boldsymbol{v}_{\\rm mesh}` is the mesh velocity
        * :math:`\\Delta t` is the time step length

        If ``tidal_dissipation = True``, the corresponding volumetric heating term defined in :func:`m_rheology.tidal_heating` is added

        .. math::
        
            - \\int_\\Omega H_{\\rm tidal}(\\eta)\\ \\hat{T}\\ \\textrm{d}\\Omega

        In case of radiation boundary condition at the top boundary ``BC_T_top = "radiation"``, the corresponding boundary condition is added

        .. math::
        
            - \\int_{\\Gamma_{\\rm top}} (1 - A)I\\hat{T}\\ + T^4\\varepsilon\\sigma_{SB}\\hat{T}\\ \\textrm{d}\\Gamma_{\\rm top}

        where

        * :math:`I` is the insolation at the body's surface (in :mod:`m_parameters`\ )
        * :math:`A` is the albedo of the body's surface (in :mod:`m_parameters`\ )
        * :math:`\\varepsilon` is the body's surface emmisivity (in :mod:`m_parameters`\ )
        * :math:`\\sigma_{SB}` is the Stefan-Boltzmann constant (in :mod:`m_constants`\ )
        
        .. note::

            If ``BC_heat_transfer[0][0] = "radiation"``, nonlinear solver will be automatically used.

        """
        if (nonlinear_heat_equation == True):
            eq_energy =  (rho_s*cp(self.Temp, self.composition)*(self.Temp - self.Temp_k)/self.dt*self.Temp_
                        
                    + 0.5*(rho_s*cp(self.Temp, self.composition)*dot(self.v_k - self.v_mesh, nabla_grad(self.Temp))*self.Temp_ 
                        + k(self.Temp_k, self.composition)*dot(nabla_grad(self.Temp), nabla_grad(self.Temp_)))

                    + 0.5*(rho_s*cp(self.Temp_k, self.composition)*dot(self.v_k - self.v_mesh, nabla_grad(self.Temp_k))*self.Temp_ 
                        + k(self.Temp_k, self.composition)*dot(nabla_grad(self.Temp_k), nabla_grad(self.Temp_))))*dx

            if (tidal_dissipation == True):
                    eq_energy += - tidal_heating(self.visc)*self.Temp_*dx

            if (BC_heat_transfer[0][0] == "radiation"):
                eq_energy += - (1.0 - albedo)*insolation*self.Temp_*self.ds(1)\
                            + pow(self.Temp, 4)*emis*SB_const*self.Temp_*self.ds(1)

            J_eq_energy = derivative(eq_energy, self.Temp)
            problem_energy = NonlinearVariationalProblem(eq_energy, self.Temp, self.bc_temp, J_eq_energy)
            self.solver_energy  = NonlinearVariationalSolver(problem_energy)

            prm = self.solver_energy.parameters
            prm["newton_solver"]["absolute_tolerance"] = 1E-8
            prm["newton_solver"]["relative_tolerance"] = 1E-7
            prm["newton_solver"]["maximum_iterations"] = 25
            prm["newton_solver"]["relaxation_parameter"] = 1.0
            prm["newton_solver"]["linear_solver"] = "bicgstab"
        
        else:

            eq_energy = (rho_s*cp(self.Temp_k, self.composition)*(self._Temp - self.Temp_k)/self.dt*self.Temp_
                        
                    + 0.5*(rho_s*cp(self.Temp_k, self.composition)*dot(self.v_k - self.v_mesh, nabla_grad(self._Temp))*self.Temp_ 
                        + k(self.Temp_k, self.composition)*dot(nabla_grad(self._Temp), nabla_grad(self.Temp_)))

                    + 0.5*(rho_s*cp(self.Temp_k, self.composition)*dot(self.v_k - self.v_mesh, nabla_grad(self.Temp_k))*self.Temp_ 
                        + k(self.Temp_k, self.composition)*dot(nabla_grad(self.Temp_k), nabla_grad(self.Temp_))))*dx
            
            if (tidal_dissipation == True):
                    eq_energy += -tidal_heating(self.visc)*self.Temp_*dx

            problem_energy = LinearVariationalProblem(lhs(eq_energy), rhs(eq_energy), self.Temp, self.bc_temp)
            self.solver_energy  = LinearVariationalSolver(problem_energy)

    def equation_Stokes(self):
        """ 
        :var: Defines the weak form of the Stokes problem

        :returns:

        Weak form of the continuity equation 

        .. math::
        
            0 = \\int_\\Omega \\nabla\\cdot \\boldsymbol{v}\\ \\hat{p}\\ \\textrm{d}\\Omega    

        Weak form of the momentum equation if ``elasticity = False``

        .. math::
        
            0 = \\int_\\Omega \\hat{p}\\nabla\\cdot \\hat{\\boldsymbol{v}} - \\rho g (\\boldsymbol{e}_z \\cdot \\hat{\\boldsymbol{v}})
            -\\eta\\ (\\nabla \\boldsymbol{v} + \\nabla^T \\boldsymbol{v}): \\nabla^T \\hat{\\boldsymbol{v}}\\ \\textrm{d}\\Omega\\\\    

        Weak form of the momentum equation if ``elasticity = True``

        .. math::
        
            0 = \\int_\\Omega \\hat{p}\\nabla\\cdot \\hat{\\boldsymbol{v}} - \\rho g (\\boldsymbol{e}_z \\cdot \\hat{\\boldsymbol{v}})
            - Z\\eta\\ (\\nabla \\boldsymbol{v} + \\nabla^T \\boldsymbol{v}): \\nabla^T \\hat{\\boldsymbol{v}}
            - (1-Z) (\\boldsymbol{\\sigma} : \\nabla^T \\hat{\\boldsymbol{v}})\\ \\textrm{d}\\Omega\\\\    

        In case of a free surface condition at the top boundary ``BC_Stokes_problem[0][0] = "free_surface"``, the 'drunken sailor' stabilization term
        is added 

        .. math::

            - \\int_{\\Gamma_{\\rm top}} \\rho_s g\\lambda\\Delta t (\\boldsymbol{v} \\cdot \\boldsymbol{n})(\\boldsymbol{e}_z \\cdot \\hat{\\boldsymbol{v}}) \\ \\textrm{d}\\Gamma_{\\rm top}

        In case of a free surface condition at the bottom boundary ``BC_Stokes_problem[1][0] = "free_surface"``, the hydrostatic pressure and the 'drunken sailor' stabilization term
        are added 

        .. math::

            + \\int_{\\Gamma_{\\rm bot}} (h_{\\rm bot}\\rho_l - h\\rho_s)(\\boldsymbol{n} \\cdot \\hat{\\boldsymbol{v}}) + (\\rho_l - \\rho_s) g\\lambda\\Delta t [(\\boldsymbol{v} + \\boldsymbol{u})
            \\cdot \\boldsymbol{n}](\\boldsymbol{e}_z \\cdot \\hat{\\boldsymbol{v}}) \\ \\textrm{d}\\Gamma_{\\rm bot}

        """
        eq_cont = div(self._v)*self.p_*dx

        if (elasticity == False):
            eq_momentum = (self._p*div(self.v_)
                - rho(self.Temp, self.composition, self.xm)*g*dot(self.e_z, self.v_)
                - self.visc*inner(2.0*sym(nabla_grad(self._v)), nabla_grad(self.v_).T))*dx     
            
        else:

            eq_momentum = (self._p*div(self.v_)
                - rho(self.Temp, self.composition, self.xm)*g*dot(self.e_z, self.v_)
                - z(self.visc, self.shear_modulus, self.dt)*self.visc*inner(2.0*sym(nabla_grad(self._v)), nabla_grad(self.v_).T)\
                - (1.0 - z(self.visc, self.shear_modulus, self.dt))*inner(self.stress_dev_tensor, nabla_grad(self.v_).T))*dx
        
        # --- Top surface "drunken sailor" stabilization ---
        if (BC_Stokes_problem[0][0] == "free_surface"): 
            eq_momentum += -rho_s*g*self._lambda*self.dt*dot(self._v, self.normal)*dot(self.e_z, self.v_)*self.ds(1)
            
        # --- Hydrostatic pressure and bottom surface "drunken sailor"---
        if (BC_Stokes_problem[1][0] == "free_surface" or phase_transition == True): 
            eq_momentum +=  dot(self.normal, self.v_)*(self.h_bot*rho_l - height*rho_s)*g*self.ds(2)\
            + (rho_l - rho_s)*g*self._lambda*self.dt*dot(self._v + self.u, self.normal)*dot(self.e_z, self.v_)*self.ds(2)
             
        problem_stokes = LinearVariationalProblem(lhs(eq_cont + eq_momentum),rhs(eq_cont + eq_momentum), self.pv, self.bc_stokes)
        self.solver_stokes = LinearVariationalSolver(problem_stokes)

    def equation_free_surface_top(self):
        """ 
        :var: Defines the weak form, ``LinearVariationalProblem`` and ``LinearVariationalSolver`` for the top surface evolution equation.
        :returns: 

        .. math::
            \\begin{align}
            0 = &\\int_\\Omega \\frac{\\partial h_{\\rm top}}{\\partial z} \\frac{\\partial \\hat{h}_{\\rm top}}{\partial z}\\ \\textrm{d}\\Omega\\\\
            &- \\int_{\\Gamma_{\\rm top}} \\frac{\\partial h_{\\rm top}}{\\partial z} \\hat{h}_{\\rm top}\\ \\textrm{d}\\Gamma_{\\rm top}\\\\
            &- \\int_{\\Gamma_{\\rm top}} \\left[h_{\\rm top} - h_{\\rm top}^k - \\Delta t \\left(\\frac{\\partial h_{\\rm top}}{\\partial x} v_x - v_z \\right)\\right]
            \\left(\\frac{\\partial \\hat{h}_{\\rm top}}{\\partial z} - \\frac{\\hat{h}_{\\rm top}}{\\gamma h_e} \\right)\\ \\textrm{d}\\Gamma_{\\rm top}
            \\end{align}

        where

        * :math:`h_{\\rm top}` is the top topography trial function
        * :math:`h_{\\rm top}^k` is the top topography function from the previous time step
        * :math:`\\hat{h}_{\\rm top}` is the top topography test function
        * :math:`v_x` and :math:`v_z` are the horizontal and vertical components of the velocity
        * :math:`\\Delta t` is the time step length
        * :math:`\\gamma` is a stabilization parameter
        * :math:`h_e` is the element size at the top boundary

        """

        # --- Top free surface equation ---
        eq_fs_top = (self._h_top.dx(1)*self.h_top_.dx(1))*dx\
            - (self._h_top.dx(1)*self.h_top_)*self.ds(1)\
            - ((self._h_top - self.h_top_k + self.dt*(self._h_top.dx(0)*self.v_k[0] - self.v_k[1]))\
            * (self.h_top_.dx(1) - self.h_top_/(self.gamma*self.h_e_top)))*self.ds(1)

        bc_fs_top = []

        free_surface_top = LinearVariationalProblem(lhs(eq_fs_top), rhs(eq_fs_top), self.h_top, bc_fs_top)
        self.solver_free_surface_top = LinearVariationalSolver(free_surface_top)

    def equation_free_surface_bottom(self):
        """ 
        :var: Defines the weak form, ``LinearVariationalProblem`` and ``LinearVariationalSolver`` for the bottom surface evolution equation.
        :returns: 

        .. math::
            \\begin{align}
            0 = &\\int_\\Omega \\frac{\\partial h_{\\rm bot}}{\\partial z} \\frac{\\partial \\hat{h}_{\\rm bot}}{\partial z}\\ \\textrm{d}\\Omega\\\\
            &- \\int_{\\Gamma_{\\rm bot}} \\frac{\\partial h_{\\rm bot}}{\\partial z} \\hat{h}_{\\rm bot}\\ \\textrm{d}\\Gamma_{\\rm bot}\\\\
            &- \\int_{\\Gamma_{\\rm bot}} \\left[h_{\\rm bot} - h_{\\rm bot}^k - \\Delta t \\left(\\frac{\\partial h_{\\rm bot}}{\\partial x} (v_x + u_x) - (v_z + u_z) \\right)\\right]
            \\left(\\frac{\\partial \\hat{h}_{\\rm bot}}{\\partial z} - \\frac{\\hat{h}_{\\rm bot}}{\\gamma h_e} \\right)\\ \\textrm{d}\\Gamma_{\\rm bot}
            \\end{align}

        where

        * :math:`h_{\\rm bot}` is the bottom topography trial function
        * :math:`h_{\\rm bot}^k` is the bottom topography function from the previous time step
        * :math:`\\hat{h}_{\\rm bot}` is the bottom topography test function
        * :math:`v_x` and :math:`v_z` are the horizontal and vertical components of the velocity
        * :math:`u_x` and :math:`u_z` are the horizontal and vertical components of the phase change velocity
        * :math:`\\Delta t` is the time step length
        * :math:`\\gamma` is a stabilization parameter
        * :math:`h_e` is the element size at the bottom boundary

        """

        # --- Bottom free surface equation ---
        eq_fs_bot = (self._h_bot.dx(1))*(self.h_bot_.dx(1))*dx\
            - (self._h_bot.dx(1))*self.h_bot_*self.ds(2)\
            - (self._h_bot - self.h_bot_k + self.dt*(self._h_bot.dx(0)*(self.v_k[0] + self.u[0]) - (self.v_k[1] + self.u[1])))\
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

    def equation_mesh_displacement(self):
        """ 
        :var: Defines the weak form, ``LinearVariationalProblem`` and ``LinearVariationalSolver`` for mesh displacement distribution equation.
        :returns: 

        If ``mesh_displacement_laplace = "full"`` then

        .. math::
        
            0 = \\int_\\Omega \\nabla h_2 \\cdot \\nabla \\hat{h}_2\\ \\textrm{d}\\Omega

        If ``mesh_displacement_laplace = "z-only"`` then

        .. math::
        
            0 = \\int_\\Omega \\frac{\\partial h_2}{\\partial z}\\frac{\\partial \\hat{h}_2}{\\partial z} \\ \\textrm{d}\\Omega
        
        where

        * :math:`h_2` is the mesh displacement trial function
        * :math:`\\hat{h}_2` is the mesh displacement test function

        """

        # --- Harmonics mesh displacement distribution ---
        if (mesh_displacement_laplace == "full"):
            eq_surf_move = dot(nabla_grad(self._h2), nabla_grad(self.h2_))*dx

        if (mesh_displacement_laplace == "z-only"):
            eq_surf_move = (self._h2.dx(1) * self.h2_.dx(1))*dx 

        bc_sm = [DirichletBC(self.sCG1, self.h2_top, self.boundary_parts, 1),\
                 DirichletBC(self.sCG1, self.h2_bot, self.boundary_parts, 2)]

        problem_surf_move = LinearVariationalProblem(lhs(eq_surf_move), rhs(eq_surf_move), self.h2, bc_sm)
        self.solver_surf_move = LinearVariationalSolver(problem_surf_move)

    def solve_initial_heat_equation(self):
        
        if (initial_tidal_dissipation == True):
            try:
                self.visc.assign(project(eta_eff(self.Temp, self.xm), self.sDG0))
                self.solver_energy_ini.solve()
                
            except:
                if (rank == 0):
                    print("\n\tError:")
                    print("\tPrescribed amplitude of tidal heating might not compatible with conductive state.\n")
                exit()

        else:
            self.solver_energy_ini.solve()

        self.q_cond_top.assign(assemble(dot(-k(self.Temp, self.composition)*nabla_grad(self.Temp), self.normal)*self.ds(1)) / self.top_length)
        self.q_ice_aver_ini.assign(assemble(sqrt(dot(self.q_ice, self.q_ice))*self.ds(2))/assemble(self.unit_scalar*self.ds(2)))
        self.Temp_k.assign(self.Temp)
    
    def solve_heat_equation(self):
        if (internal_melting == True):
            # --- Algorithm from Kalousova (2015) ---
            if (self.only_melt_tracers == True):
                self.add_melt_tracers()

            self.solver_energy.solve()

            melt_interpolation(self.mesh, self.tracers_in_cells, self.tracers, 0, self.xm_k)

            # Compute change in porosity and in temperature
            self.compute_melting()

            # Update temperature and porosity
            melt_interpolation(self.mesh, self.tracers_in_cells, self.tracers, 1, self.delta_T)
            self.Temp.assign(project(self.Temp + self.delta_T, self.sCG2))
            
            melt_interpolation(self.mesh, self.tracers_in_cells, self.tracers, 0, self.xm)
            
            # Evaluate melting/crystallization rate
            self.melting_rate.assign(project(rho_m*(self.xm - self.xm_k)/self.dt*rho_s/rho_m, self.sDG0))

            if (self.only_melt_tracers == True):
                self.delete_melt_tracers()            
            # -----------------
        else:
            self.solver_energy.solve()

        self.Temp_k.assign(self.Temp)
        self.q_ice.assign(project(-k(self.Temp, self.composition)*nabla_grad(self.Temp), self.vCG1))

    def solve_Stokes_problem(self, step, FilesClass):
        """ Solves the Stokes problem.
        
        """
        # self.heating.assign(project(tidal_heating(self.visc), self.sDG0))

        # --- If there is the phase transition, estimate first phase change velocity
        #     as it is in the stabilization term ---
        if (BC_Stokes_problem[1][0] == "free_surface" and phase_transition == True):
            self.compute_u()
        
        Picard_iter = 0
        v_error = 1.0

        # --- After the tracers are advected (from the previous loop) ---
        # scalar_interpolation(self.mesh, self.tracers_in_cells, self.tracers, 6, "ARITM", self.plastic_strain)

        if (elasticity == True):
            self.shear_modulus.assign(project(G(self.composition), self.sDG0))

        self.update_viscosity(step, Picard_iter)

        if (elasticity == True):
            self.z_function.assign(project(z(self.visc, self.shear_modulus, self.dt), self.sDG0))

        # --- Picard iterations of Stokes solver ---  
        while (Picard_iter < Picard_iter_max and v_error > Picard_iter_error):

            self.solver_stokes.solve()        
            if (stokes_elements == "Mini"):
                self.p_out, self.v_out, self.vb = self.pv.split(True)

            if (stokes_elements == "TH"):
                self.p_out, self.v_out = self.pv.split(True)
        

            if (BC_Stokes_problem[0][0] == "free_surface" and BC_Stokes_problem[1][0] == "free_surface"):

                if (stokes_null == "volume"):
                    self.v_aver.assign(assemble(self.v_out[1]*dx)/assemble(self.unit_scalar*dx))
                    self.v_out.assign(project(self.v_out - self.v_aver*self.e_z, self.vCG1))

                if (stokes_null == "boundary"):
                    self.v_top_aver.assign(assemble(self.v_out[1]*self.ds(1))/assemble(self.unit_scalar*self.ds(1)))
                    self.v_out.assign(project(self.v_out - self.v_top_aver*self.e_z, self.vCG1))

            if (self.nonlinear_rheology == True):
                self.iteration_error.assign(project(sqrt(dot(self.v_out - self.v_k, self.v_out - self.v_k))/sqrt(dot(self.v_k, self.v_k)), self.sCG1))
        
                if (step == 0 and Picard_iter == 0):
                    v_error = 1.0
                else:
                    if (error_type == "maximum"):
                        v_error = MPI.max(comm, np.abs(self.iteration_error.vector()).max()) 

                    if (error_type == "integrated"):
                        v_error = assemble(self.iteration_error*dx)/assemble(self.unit_scalar*dx)

                if (rank == 0):
                    print("\n\tPicard iteration %d. Error = %6.3E\n" % (Picard_iter, v_error))

                    file = open("data_" + self.name + "/picard_iter.dat","a")
                    file.write(("%.5E\t"+"%.5E"+"\n")%(step + Picard_iter/Picard_iter_max, v_error))    
                    file.close() 
                
                self.strain_rate_inv.assign(project(strain_rate_II(self.v_out), self.sDG0))
                
                self.p_k.assign(self.p_out)
                self.v_k.assign(self.v_out)

                if (elasticity == True and viscosity_type == "GK_2001"):
                    get_new_stress_iter(self.mesh, self.Temp, self.xm, self.stress_dev_inv,\
                                        self.stress_dev_inv_k, self.strain_rate_inv, self.composition,\
                                        self.shear_modulus, step, Picard_iter, self.dt)

                Picard_iter += 1
                
                self.update_viscosity(step, Picard_iter)
                if (elasticity == True):
                    self.z_function.assign(project(z(self.visc, self.shear_modulus, self.dt), self.sDG0))

            else:
                self.strain_rate_inv.assign(project(strain_rate_II(self.v_out), self.sDG0))
                if (elasticity == True):
                    self.z_function.assign(project(z(self.visc, self.shear_modulus, self.dt), self.sDG0))

                self.p_k.assign(self.p_out)
                self.v_k.assign(self.v_out)
                break

        self.density.assign(project(rho(self.Temp, self.composition, self.xm), self.sDG0))
        if (time_step_position == "stokes"):
            self.dt.assign(time_step(self.mesh, self.v_k, self.v_mesh, H_max, self.composition, self.Temp, self.unit_scalar, self.t))        

        # --- Estimate the visco-elastic stress ---
        get_new_stress(self.mesh, self.Temp, self.xm, self.stress_dev_inv,\
                    self.stress_dev_inv_k, self.strain_rate_inv, self.composition,\
                    self.shear_modulus, step, Picard_iter, self.dt)

        if (plasticity == True):
            self.yield_stress.assign(project(sigma_yield(self.p_k, self.plastic_strain, self.composition), self.sDG0)) 
            self.yield_function.assign(project(self.stress_dev_inv - self.yield_stress, self.sDG0))

            # --- Integrate plastic strain on tracers ---
            plastic_strain_integration(self.mesh, self.tracers_in_cells, self.tracers, self.dt, self.yield_function, self.strain_rate_inv)
            
    def update_stress(self):
        # --- Update stress on tracers ---
        self.strain_rate_tensor.assign(project(sym(nabla_grad(self.v_k)), self.tDG0))
        stress_update(self.mesh, self.tracers_in_cells, self.tracers, self.visc, self.strain_rate_tensor, self.z_function)

        if (plasticity == True):
            # --- Reduce stress where the ice is yielding ---
            stress_reduction(self.mesh, self.tracers_in_cells, self.tracers, self.yield_function, self.yield_stress, self.stress_dev_inv)

    def rotate_and_interpolate_stress(self):
        # --- Rotate stress ---
        self.vorticity_tensor.assign(project(skew(nabla_grad(self.v_k)), self.tDG0))
        stress_rotation(self.mesh, self.tracers_in_cells, self.tracers, self.vorticity_tensor, self.dt)

        stress_interpolation(self.mesh, self.tracers_in_cells, self.tracers, self.stress_dev_tensor)
        self.stress_dev_inv_k.assign(project(tensor_2nd_invariant(self.stress_dev_tensor), self.sDG0))

        if (rank==0):
            print("Done.")

    def update_viscosity(self, step, Picard_iter):
        """ 
        :var: Updates the viscosity fuction 
        """
        self.sr_min.assign(MPI.min(comm, self.strain_rate_inv.vector().min()))

        self.visc.assign(project(eta_eff(self.p_k, self.Temp, self.strain_rate_inv,\
                                        self.stress_dev_inv, self.xm, self.composition, self.plastic_strain,\
                                        step, Picard_iter, self.stress_dev_inv_k, self.shear_modulus, self.dt, self.sr_min), self.sDG0))
        
        self.shear_modulus.assign(project(G(self.composition), self.sDG0))

    def compute_u(self):
        """ 
        :var: Updates the heat fluxes :math:`\\boldsymbol{q}_{\\rm ice}`  and :math:`\\boldsymbol{q}_{\\rm ocean}` and computes
              new phase transition velocity :math:`\\boldsymbol{u}`\ .
        """
        self.q_ice.assign(project(-k(self.Temp, self.composition)*nabla_grad(self.Temp), self.vCG1))
        self.q_ice_aver.assign(assemble(sqrt(dot(self.q_ice, self.q_ice))*self.ds(2))/assemble(self.unit_scalar*self.ds(2)))
        self.q_water.assign(project((self.q_ice_aver - DAL_factor*self.h_bot)*self.e_z, self.vCG1))
        self.u.assign(project(- dot(self.n_bot, self.q_ice - self.q_water)*self.n_bot/(Lt*rho_s), self.vCG1))

    def solve_topography_evolution(self):
        # --- Compute the top topography evolution ---
        if (BC_Stokes_problem[0][0] == "free_surface"):
            self.solver_free_surface_top.solve()

        # --- Upwards-facing normal vector to the bottom boundary ---
        self.n_bot.assign(project((self.e_z - self.e_x*dot(nabla_grad(self.h_bot), self.e_x))/sqrt(1.0 + dot(nabla_grad(self.h_bot), self.e_x)**2), self.vCG1))

        # --- Compute the velocity of the phase change ---
        if (BC_Stokes_problem[1][0] == "free_surface" and phase_transition == True):
            self.compute_u()

        # --- Compute the bottom topography evolution ---
        if (BC_Stokes_problem[1][0] == "free_surface"):
            self.solver_free_surface_bot.solve()

        # --- Compute the mesh displecement distribution ---
        if (BC_Stokes_problem[0][0] == "free_surface" or BC_Stokes_problem[1][0] == "free_surface"):    
            self.h2_top.assign(self.h_top - self.h_top_k)
            self.h2_bot.assign(self.h_bot - self.h_bot_k)

            self.solver_surf_move.solve()
            self.u_mesh.assign(project(self.h2*self.e_z, self.vCG1))

            self.h_top_k.assign(self.h_top)
            self.h_bot_k.assign(self.h_bot)

            self.v_mesh.assign(project(self.u_mesh/self.dt, self.vCG1))

        if (time_step_position == "end"):
            self.dt.assign(time_step(self.mesh, self.v_k, self.v_mesh, H_max, self.composition, self.Temp, self.unit_scalar, self.t))    