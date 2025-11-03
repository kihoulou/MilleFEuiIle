from dolfin import *
import numpy
from scipy import special
from m_parameters import *

def tidal_heating(visc):

    # --- Set for Europa ---
    visc_max = 1.5e14
    nu = visc_max/visc

    if (heating_model == "Maxwell"):
        H_value = 2.0*H_max/(nu + 1.0/nu)

    elif (heating_model == "Andrade"):
        gamma_tidal = special.gamma(1.0 + alpha_and)

        H_value = 2.0*H_max*(1.0 + nu**(alpha_and-1.0)*gamma_tidal*sin(alpha_and*pi/2.0))\
        /(nu + 1.0/nu + nu**(alpha_and-1.0)*gamma_tidal*(nu**alpha_and*gamma_tidal + 2.0*cos(alpha_and*pi/2.0) + 2.0*nu*sin(alpha_and*pi/2.0)))

    return H_value

def tensor_2nd_invariant(tensor):
    return sqrt(0.5*inner(tensor, tensor))

def strain_rate_II(v):

    return sqrt(inner(sym(nabla_grad(v)), sym(nabla_grad(v))) / 2.0)

def max_function(a, b):

    return (a + b)/2.0 + abs(a - b)/2.0

# --- Plasticity ---
def cohesion(plastic_strain):
    value = cohesion_strong + (cohesion_weak - cohesion_strong)*plastic_strain/eps_weak

    return conditional(lt(plastic_strain, eps_weak), value, cohesion_weak)

def sigma_yield(p_k, plastic_strain, composition):

    angle = int_friction_angle*(np.pi/180.0)
    value = p_k*sin(angle) + cohesion(plastic_strain)*cos(angle)

    # angle = (int_friction_angle*composition[0] + int_friction_angle2*composition[1])*(np.pi/180.0)
    # value = p_k*sin(angle) + cohesion(plastic_strain)*cos(angle)

    # friction_coef = 0.6
    # value = (cohesion(plastic_strain) + friction_coef*(p_k))*cos(atan(friction_coef))

    return conditional(lt(value, yield_stress_min), yield_stress_min, conditional(gt(value, yield_stress_max), yield_stress_max, value))

def plastic_strain_integration(mesh, tracers_in_cells, tracers, dt, yield_function, strain_rate):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]
            px = tracers[tracer_no][0]
            py = tracers[tracer_no][1]

            yield_function_tracer = yield_function(Point(px,py))
            strain_rate_tracer = strain_rate(Point(px,py))

            if (yield_function_tracer >= 0.0):
                tracers[tracer_no][6] += float(dt)*strain_rate_tracer

            if (healing == True):
                tracers[tracer_no][6] = tracers[tracer_no][6]/(1.0 + float(dt)/healing_timescale)

# --- Elasticity ---
def z(visc, shear_modulus, dt):
    return dt/(dt + visc/shear_modulus)

def get_new_stress(mesh, Temp, xm, stress_dev_inv, stress_dev_inv_k, strain_rate_invariant, composition, shear_modulus, step, Picard_iter, dt):
    ranks = []
    eval_type = "local"
    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()
        
        temp_cell = Temp(Point(centroid.x(), centroid.y()))
        sigma_II_cell = stress_dev_inv(Point(centroid.x(), centroid.y()))
        eps_II_cell = strain_rate_invariant(Point(centroid.x(), centroid.y()))
        
        
        xm_cell = xm(Point(centroid.x(), centroid.y()))
        comp_cell = []
        for i in range (len(composition)):
            comp_cell.append(composition[i](Point(centroid.x(), centroid.y())))
        
        visc_cell = eta_ductile(temp_cell, eps_II_cell, sigma_II_cell, xm_cell, comp_cell, step, Picard_iter, eval_type)
        
        if (elasticity == True):
            sigma_II_k_cell = stress_dev_inv_k(Point(centroid.x(), centroid.y()))
            shear_mod_cell = shear_modulus(Point(centroid.x(), centroid.y()))
            z_cell = z(visc_cell, shear_mod_cell, dt)

            new_stress = 2.0*z_cell*visc_cell*eps_II_cell + (1.0 - z_cell)*sigma_II_k_cell
        else:
            new_stress = 2.0*visc_cell*eps_II_cell

        ranks.append(new_stress)
        
    ranks = numpy.array(ranks)                                                            
    stress_dev_inv.vector().set_local(ranks)

def get_new_stress_iter(mesh, Temp, xm, stress_dev_inv, stress_dev_inv_k, strain_rate_inv, composition, shear_modulus, step, Picard_iter, dt):
    ranks = []
    eval_type = "local"

    for j in range(mesh.num_cells()):
        centroid    = Cell(mesh, j).midpoint()

        temp_cell       = Temp(Point(centroid.x(), centroid.y()))
        sigma_II_cell   = stress_dev_inv(Point(centroid.x(), centroid.y()))
        sigma_II_k_cell = stress_dev_inv_k(Point(centroid.x(), centroid.y()))
        eps_II_cell     = strain_rate_inv(Point(centroid.x(), centroid.y()))

        xm_cell = xm(Point(centroid.x(),centroid.y()))
        sm_cell = shear_modulus(Point(centroid.x(),centroid.y()))

        comp_cell = []
        for i in range (len(composition)):
            comp_cell.append(composition[i](Point(centroid.x(), centroid.y())))

        visc_cell = eta_ductile(temp_cell, eps_II_cell, sigma_II_cell, xm_cell, comp_cell, step, Picard_iter, eval_type)

        stress_l = sigma_II_k_cell
        stress_r = sigma_II_k_cell
        
        z_cell = float(dt)/(float(dt) + visc_cell/sm_cell)
        
        new_stress = 2.0*z_cell*visc_cell*eps_II_cell + (1.0 - z_cell)*sigma_II_k_cell
        
        if (new_stress < sigma_II_k_cell):
            stress_l = new_stress
        if (new_stress > sigma_II_k_cell):
            stress_r = new_stress

        visc_cell = eta_ductile(temp_cell, eps_II_cell, new_stress, xm_cell, comp_cell, step, Picard_iter, eval_type)

        z_cell = float(dt)/(float(dt) + visc_cell/sm_cell)
        new_stress = 2.0*z_cell*visc_cell*eps_II_cell + (1.0-z_cell)*sigma_II_k_cell

        if (new_stress < stress_l):
            stress_l = new_stress
        if (new_stress > stress_r):
            stress_r = new_stress

        # --- We have left and right boundary, let's iterate ---
        error = 1
        k = 0

        while (error > stress_iter_error and k < 50):
            stress_c = (stress_l + stress_r)/2.0
            visc_cell = eta_ductile(temp_cell, eps_II_cell, stress_c, xm_cell, comp_cell, step, Picard_iter, eval_type)

            z_cell = float(dt)/(float(dt) + visc_cell/sm_cell)
            new_stress = 2.0*z_cell*visc_cell*eps_II_cell + (1.0 - z_cell)*sigma_II_k_cell

            # Criteria based only on the invariant
            if (new_stress < stress_c):
                stress_r = stress_c

            if (new_stress > stress_c):
                stress_l = stress_c

            error = abs((new_stress-stress_c)/stress_c)
            k += 1

        ranks.append(stress_c)

    MPI.barrier(comm)
    ranks = numpy.array(ranks)
    stress_dev_inv.vector().set_local(ranks)

# def get_mechanisms(mesh,Temp,stress_invariant,old_stress_invariant,strain_rate_invariant,mechanisms,p_k, yield_stress, dt, sr_min):
#     ranks = []
#     eval_type = "local"

#     for j in range(mesh.num_cells()):
#         centroid = Cell(mesh, j).midpoint()

#         temp_cell = Temp(Point(centroid.x(),centroid.y()))
#         eps_cell = strain_rate_invariant(Point(centroid.x(),centroid.y()))
#         tau_cell = stress_invariant(Point(centroid.x(),centroid.y()))
#         tau_old_cell = old_stress_invariant(Point(centroid.x(),centroid.y()))
#         p_cell = p_k(Point(centroid.x(),centroid.y()))
#         ys_cell = yield_stress(Point(centroid.x(),centroid.y()))

#         visc_diff = eta_diff(temp_cell)
#         visc_disl = eta_disl(temp_cell, eps_cell, tau_cell, eval_type)
#         visc_GBS_BS = eta_GBS(temp_cell, eps_cell, tau_cell,eval_type) + eta_BS(temp_cell, eps_cell, tau_cell, p_cell)
        
#         if (plasticity==True):
#             visc_plast = 0.5*ys_cell/max_function(eps_cell - (ys_cell - tau_old_cell)/(2.0*shear_modulus*float(dt)), float(sr_min)*1e-6)
#             array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, visc_plast],[5.0, eta_max]]
#         else:
#             array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, eta_max]]

#         array_mech.sort(key=lambda x: x[1])

#         ranks.append(array_mech[0][0])
    
#     ranks = numpy.array(ranks)
#     mechanisms.vector().set_local(ranks)

def G(composition):

    # --- Rising plume benchmark ---
    # return G_mantle**composition[0] * G_lid**composition[1] * G_plume**composition[2]

    # --- Ice shell ---
    return 3.52e9

def eta_ductile(Temp, strain_rate, stress, xm, composition, step, Picard_iter, eval_type):
    
    if (viscosity_type=="constant"):
        eta_v = eta_0

    if (viscosity_type == "temp-dep"):
        T_bot = BC_heat_transfer[1][1]
        eta_v = 1.0/(1.0/(eta_0*exp(Q_activ/R_gas*(1.0/Temp - 1.0/T_bot))) + 1.0/eta_max)*exp(-45.0*xm)
        # eta_v = eta_0*exp(-14.0*(Temp - T_bot)/170.0)*exp(-45.0*xm)

    if (viscosity_type == "GK_2001"):
        if (step == 0 and Picard_iter == 0):
            eta_v =  eta_diff(Temp)*exp(-45.0*xm) 

        else:               
            eta_v = 1.0/(1.0/eta_diff(Temp)\
                    + 1.0/eta_disl(Temp, strain_rate, stress, eval_type)\
                    + 1.0/(eta_GBS(Temp, strain_rate, stress, eval_type) + eta_BS(Temp, strain_rate, stress))\
                    )*exp(-45.0*xm)
            
    if (viscosity_type == "composition"):
        # --- Shear bands benchmark ---
        # eta_v = 1e5**composition[0]*1.0**composition[1]

        # --- Rising plume benchmark ---
        # eta_v = eta_mantle**composition[0] * eta_lid**composition[1] * eta_plume**composition[2]
        pass

    return eta_v

def eta_eff(p_k, Temp, strain_rate, stress, xm, composition, plastic_strain, step, Picard_iter, stress_dev_inv_k, shear_modulus, dt, sr_min):
    
    eval_type = "mesh"

    # --- Ductile viscosity ---
    eta_v = eta_ductile(Temp, strain_rate, stress, xm, composition, step, Picard_iter, eval_type)

    if (plasticity == True):

        if (step == 0 and Picard_iter == 0):
            # --- Visco-(elasto)-plastic rheology in the very first time step ---
            return conditional(lt(eta_v, eta_max), eta_v, eta_max)

        else:
            
            # --- Plastic viscosity ---
            if (elasticity == False or (elasticity == True and step == 0)):
                eta_p = 0.5*sigma_yield(p_k, plastic_strain, composition)/strain_rate
                
            else:
                eta_p = 0.5*sigma_yield(p_k, plastic_strain, composition)/max_function(strain_rate\
                - (sigma_yield(p_k, plastic_strain, composition) - stress_dev_inv_k)/(2.0*shear_modulus*dt), sr_min*1e-6)

            # --- Visco-(elasto)-plastic rheology after the very first time step ---
            return conditional(lt(eta_p, eta_min_plast), eta_min_plast, conditional(lt(eta_p, eta_v), eta_p, eta_v))

            # This creates sharper faults ?
            # return 1.0/(1.0/eta_v + 1.0/(eta_p + eta_min_plast) + 1.0/eta_max)

    else:
        # --- Viscous rheology only ---
        if (lower_cutoff == False and upper_cutoff == False):
            return eta_v
        
        elif (lower_cutoff == True and upper_cutoff == False):
            return conditional(lt(eta_v, eta_min), eta_min, eta_v)

        elif (lower_cutoff == False and upper_cutoff == True):
            return conditional(lt(eta_v, eta_max), eta_v, eta_max)
        
        else:
            return conditional(lt(eta_v, eta_max), conditional(lt(eta_v, eta_min), eta_min, eta_v), eta_max)
        
def eta(Q, A, n, m, Temp, strain_rate_inv, stress_inv):

    # See Gerya (2009) chapter 6.2 for derivation of prefactors
    # The correction applies also to diffusion creep, see the function below and eq. (6.8b) in Gerya (2009)

    # --- Strain rate-based formula ---
    if (elasticity == False):
        return 1.0/(3.0**((n+1)/(2.0*n))*2.0**((n-1.0)/n))*A**(-1.0/n)*d_grain**(m/n)*strain_rate_inv**((1.0-n)/n)*exp(Q/(n*R_gas*Temp))
        
    # --- Stress-based formula ---
    if (elasticity == True):
        return 1.0/(3.0**((n+1)/2.0))*A**(-1.0)*d_grain**m*stress_inv**(1.0-n)*exp(Q/(R_gas*Temp))

def eta_diff(Temp):
    V_m = 1.97e-5           #m^3
    D_0v = 9.10e-4          #m^2 s^-1
    D_0b = 6.4e-4           #m^2 s^-1
    Q_v = 59.4e3            #J mol^-1
    Q_b = 49.0e3            #J mol^-1
    delta= 9.04e-10         #m

    D_v = D_0v*exp(-Q_v/(R_gas*Temp))
    D_b = D_0b*exp(-Q_b/(R_gas*Temp))

    # 3/2 come from the scaling between diff. stress, see Eq. (6.8b) in Gerya (2009)
    return 1.0/(2*(3.0/2.0*14)*V_m)*R_gas*Temp*d_grain**2/(D_v + np.pi*delta*D_b/gs)
    
def eta_disl(Temp, strain_rate_II, stress_II, eval_type):
    n = 4.0           # -
    p = 0.0           # -
    T_crit = 258      # K
    
    # Below critical temperature
    Q_below = 60e3         #J/mol
    A_below = 4e-19        #Pa^{-4} s^{-1}

    # Above critical temperature    
    Q_above = 180e3       #J/mol
    A_above = 6e4         #Pa^{-4} s^{-1}

    if (eval_type == "mesh"):
        return conditional(gt(Temp, T_crit),\
        eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II),\
        eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II))

    if (eval_type == "local"):
        if (Temp > T_crit):
            return eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II)
        else:
            return eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II)

def eta_GBS(Temp, strain_rate_II, stress_II, eval_type):
    
    n = 1.8           # -
    p = 1.4           # -
    T_crit = 255      # K

    # Below critical temperature
    Q_below = 49e3         #J/mol
    A_below = 6.2e-14      #Pa^{-1.8} m^{1.4} s^{-1}
    
    # Above critical temperature
    Q_above = 192e3        #J/mol
    A_above = 4.8e15       #Pa^{-1.8} m^{1.4} s^{-1}

    if (eval_type == "mesh"):
        return conditional(gt(Temp, T_crit),\
        eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II),\
        eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II))  

    if (eval_type == "local"):
        if (Temp > T_crit):
            return eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II)
        else:
            return eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II)    

def eta_BS(Temp, strain_rate_II, stress_II):

    Q = 60e3          #J/mol
    A = 2.2e-7        #Pa^{-2.4} s^{-1}
    n = 2.4           # -
    p = 0.0           # -

    return eta(Q, A, n, p, Temp, strain_rate_II, stress_II)
