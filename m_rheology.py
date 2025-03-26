from dolfin import *
import numpy
from scipy import special
from m_parameters import *

def tidal_heating(visc):
    """Computes the tidal heating rate :math:`H`\ .

    :param visc: viscosity (:math:`\\eta`\ )

    :returns: 

    Defines :math:`\\nu = \\eta_{peak}/\\eta`\ , where :math:`\\eta_{peak}` is the viscosity at which the maximum of the tidal heating occurs.

    * If ``heating_model = "Maxwell"`` 
       .. math:: H = \\frac{2 H_{max}}{\\nu + 1/\\nu}

    * If ``heating_model = "Andrade"``
       .. math:: H = \\frac{2 H_{max}}{\\nu + 1/\\nu}

    """

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
    """Computes the second invariant of the strain rate tensor. 

    :param v: velocity

    :returns:
    .. math::
       \\dot{\\varepsilon}_{II} = \\sqrt{\\frac{ \\dot{\\varepsilon} : \\dot{\\varepsilon}}{2}}

    """

    return sqrt(inner(sym(nabla_grad(v)), sym(nabla_grad(v))) / 2.0)

def max_function(a,b):
    return (a+b)/2.0+abs(a-b)/2.0

# --- Plasticity ---
def cohesion(plastic_strain):
    value = cohesion_strong + (cohesion_weak - cohesion_strong)*plastic_strain/eps_weak

    return conditional(lt(plastic_strain, eps_weak), value, cohesion_weak)

def sigma_yield(p_k, plastic_strain, composition):
    angle = int_friction_angle*composition[0] + int_friction_angle2*composition[1]
    value = p_k*sin(angle) + cohesion(plastic_strain)*cos(angle)

    return conditional(lt(value, yield_stress_min), yield_stress_min, conditional(gt(value, yield_stress_max), yield_stress_max, value))

def plastic_strain_integration(mesh, tracers_in_cells, tracers, dt, yield_function, strain_rate):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]
            px = tracers[tracer_no][0]
            py = tracers[tracer_no][1]

            yield_function_tracer = yield_function(Point(px,py))
            strain_rate_tracer = strain_rate(Point(px,py))

            tol = cohesion_strong/1e9
            if (yield_function_tracer >= -tol):
                tracers[tracer_no][6] += float(dt)*strain_rate_tracer

            if (healing == True):
                tracers[tracer_no][6] = tracers[tracer_no][6]/(1.0 + float(dt)/recovery_time)

# --- Elasticity ---
def z(visc, dt):
    """Computes the viscoelastic parameter :math:`Z`\ . 

    :param visc: viscosity (:math:`\\eta`\ )
    :param dt: time step length (:math:`\\Delta t`\ )

    :returns: .. math:: Z = \\frac{\\Delta t}{\\Delta t + \\frac{\\eta}{G}}

    """
    return dt/(dt + visc/shear_modulus)

def get_new_stress_iter(mesh, Temp, xm, stress_dev_inv, stress_dev_inv_k, strain_rate_inv, step, iter, dt):
    ranks = []
    eval_type = "local"

    for j in range(mesh.num_cells()):
        centroid    = Cell(mesh, j).midpoint()
        temp_cell   = Temp(Point(centroid.x(),centroid.y()))
        
        stress_dev_inv_cell    = stress_dev_inv(Point(centroid.x(),centroid.y()))
        stress_dev_inv_k_cell  = stress_dev_inv_k(Point(centroid.x(),centroid.y()))
        strain_rate_inv_cell    = strain_rate_inv(Point(centroid.x(),centroid.y()))

        xm_cell     = xm(Point(centroid.x(),centroid.y()))

        # --- First estimate of viscosity from previous stress ---
        if (step == 0 or (step == 1 and iter == 0) and (reload_HDF5 == False or reload_tracers == False)):
            visc_cell = 1.0/(1.0/eta_diff(temp_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
        else:
            visc_cell = 1.0/(1.0/eta_diff(temp_cell)\
                            + 1.0/eta_disl(temp_cell, strain_rate_inv_cell, stress_dev_inv_cell, eval_type)\
                            + 1.0/(eta_BS(temp_cell, strain_rate_inv_cell, stress_dev_inv_cell)\
                                + eta_GBS(temp_cell, strain_rate_inv_cell, stress_dev_inv_cell, eval_type))\
                            + 1.0/eta_max)*exp(-45.0*xm_cell)

        stress_l = stress_dev_inv_k_cell
        stress_r = stress_dev_inv_k_cell

        z_cell = float(dt)/(float(dt) + visc_cell/shear_modulus)
        new_stress = 2.0*z_cell*visc_cell*strain_rate_inv_cell + (1.0-z_cell)*stress_dev_inv_k_cell

        if (new_stress < stress_dev_inv_cell):
            stress_l = new_stress
        if (new_stress > stress_dev_inv_cell):
            stress_r = new_stress

        # --- Second estimate of viscosity from the new stress ---
        if (step == 0 or (step == 1 and iter == 0)):
            visc_cell = 1.0/(1.0/eta_diff(temp_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
        else:
            visc_cell = 1.0/(1.0/eta_diff(temp_cell)\
                        + 1.0/eta_disl(temp_cell, strain_rate_inv_cell, new_stress, eval_type)\
                        + 1.0/(eta_BS(temp_cell, strain_rate_inv_cell, new_stress,)\
                               + eta_GBS(temp_cell, strain_rate_inv_cell, new_stress, eval_type))\
                         + 1.0/eta_max)*exp(-45.0*xm_cell)

        z_cell = float(dt)/(float(dt) + visc_cell/shear_modulus)
        new_stress = 2.0*z_cell*visc_cell*strain_rate_inv_cell + (1.0-z_cell)*stress_dev_inv_k_cell

        if (new_stress < stress_l):
            stress_l = new_stress
        if (new_stress > stress_r):
            stress_r = new_stress

        # --- We have left and right boundary, let's iterate ---
        error = 1
        k = 0
        # for k in range(0,iter_max):
        while (error > stress_iter_error and k < 50):
            stress_c = (stress_l + stress_r)/2.0

            if (step == 0 or (step == 1 and iter == 0)):
                visc_cell = 1.0/(1.0/eta_diff(temp_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
            else:
                visc_cell = 1.0/(1.0/eta_diff(temp_cell)\
                            + 1.0/eta_disl(temp_cell, strain_rate_inv_cell, stress_c, eval_type)\
                            + 1.0/(eta_BS(temp_cell, strain_rate_inv_cell, stress_c)\
                                   + eta_GBS(temp_cell, strain_rate_inv_cell, stress_c, eval_type))\
                            + 1.0/eta_max)*exp(-45.0*xm_cell)

            z_cell = float(dt)/(float(dt) + visc_cell/shear_modulus)
            new_stress = 2.0*z_cell*visc_cell*strain_rate_inv_cell + (1.0-z_cell)*stress_dev_inv_k_cell

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

def get_mechanisms(mesh,Temp,stress_invariant,old_stress_invariant,strain_rate_invariant,mechanisms,p_k, yield_stress, dt, sr_min):
    ranks = []
    eval_type = "local"

    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()

        temp_cell = Temp(Point(centroid.x(),centroid.y()))
        eps_cell = strain_rate_invariant(Point(centroid.x(),centroid.y()))
        tau_cell = stress_invariant(Point(centroid.x(),centroid.y()))
        tau_old_cell = old_stress_invariant(Point(centroid.x(),centroid.y()))
        p_cell = p_k(Point(centroid.x(),centroid.y()))
        ys_cell = yield_stress(Point(centroid.x(),centroid.y()))

        visc_diff = eta_diff(temp_cell)
        visc_disl = eta_disl(temp_cell, eps_cell, tau_cell, eval_type)
        visc_GBS_BS = eta_GBS(temp_cell, eps_cell, tau_cell,eval_type) + eta_BS(temp_cell, eps_cell, tau_cell, p_cell)
        
        if (plasticity==True):
            visc_plast = 0.5*ys_cell/max_function(eps_cell - (ys_cell - tau_old_cell)/(2.0*shear_modulus*float(dt)), float(sr_min)*1e-6)
            array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, visc_plast],[5.0, eta_max]]
        else:
            array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, eta_max]]

        array_mech.sort(key=lambda x: x[1])

        ranks.append(array_mech[0][0])
    
    ranks = numpy.array(ranks)
    mechanisms.vector().set_local(ranks)

def eta_visc(composition):
    return  1e5**composition[0]*1e0**composition[1]

def eta_eff(p_k, Temp, strain_rate, stress, xm, composition, plastic_strain, step, Picard_iter, stress_dev_inv_k, dt, sr_min):
    eval_type = "mesh"

    if (plasticity == True):
            if (elasticity == False):
                eta_p = 0.5*sigma_yield(p_k, plastic_strain, composition)/strain_rate
            else:
                eta_p = 0.5*sigma_yield(p_k, plastic_strain, composition)/max_function(strain_rate\
                - (sigma_yield(p_k, plastic_strain, composition) - stress_dev_inv_k)/(2.0*shear_modulus*dt), sr_min*1e-6)

    if (viscosity_type=="constant"):
        eta_v = eta_0

    if (viscosity_type == "temp-dep"):
        eta_v = 1.0/(1.0/(eta_0*exp(Q_activ/R_gas*(1.0/Temp - 1.0/T_bot))) + 1.0/eta_max)*exp(-45.0*xm)

    if (viscosity_type == "GK_2001"):
        if (step == 0 and Picard_iter == 0):
            eta_v =  1.0/(1.0/eta_diff(Temp) + 1.0/eta_max)*exp(-45.0*xm) 

        else:            
            
            eta_v = 1.0/(1.0/eta_diff(Temp)\
                    + 1.0/eta_disl(Temp, strain_rate, stress, eval_type)\
                    + 1.0/(eta_GBS(Temp, strain_rate, stress, eval_type) + eta_BS(Temp, strain_rate, stress))\
                    + 1.0/eta_max)*exp(-45.0*xm)
            
    if (viscosity_type == "composition"):
        eta_v = 1e5**composition[0]*1.0**composition[1]
            
    return conditional(lt(eta_p, eta_v), eta_p, eta_v)

    # if (plasticity == True):
    #     if (step == 0 and Picard_iter == 0):
    #         return conditional(lt(eta_v, eta_max), eta_v, eta_max)
    #     else:
    #         return conditional(lt(eta_p, eta_min_plast), eta_min_plast, conditional(lt(eta_p, eta_v), eta_p, eta_v))
    # else:
    #     return conditional(lt(eta_v, eta_max), eta_v, eta_max)

def eta(Q, A, n, m, Temp, strain_rate_inv, stress_inv):
    # See Gerya (2009) chapter 6.2 for derivation of prefactors
    # The correction applies also to diffusion creep, see the function below and eq. (6.8b) in Gerya (2009)
    # if (viscosity_from == "strain_rate"):
    if (elasticity == False):
        return 1.0/(3.0**((n+1)/(2.0*n))*2.0**((n-1.0)/n))*A**(-1.0/n)*d_grain**(m/n)*strain_rate_inv**((1.0-n)/n)*exp(Q/(n*R_gas*Temp))
        
    # if (viscosity_from == "stress"):
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
    return 1.0/(2*(3.0/2.0*14)*V_m)*R_gas*Temp*d_grain**2/(D_v + np.pi*delta*D_b/d_grain)

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
