from dolfin import *
import numpy
from scipy import special

import os 
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

    elif (heating_model == "none"):
        H_value = 0.0

    return H_value

def strain_rate_II(v):
    return sqrt(1.0/8.0*inner(nabla_grad(v)+nabla_grad(v).T,nabla_grad(v)+nabla_grad(v).T))

def max_function(a,b):
    return (a+b)/2.0+abs(a-b)/2.0

# --- Plasticity ---
def friction_angle(plastic_strain):
    value = a0 + (a_inf - a0)*plastic_strain/eps_inf
    return conditional(lt(plastic_strain, eps_inf), value, a_inf)

def cohesion(plastic_strain):
    value = C0 + (C_inf - C0)*plastic_strain/eps_inf
    return conditional(lt(plastic_strain, eps_inf), value, C_inf)

def yield_criterion(yield_stress, stress_invariant):
    return stress_invariant - yield_stress

def angle_friction(composition):
    # --- Shear bands benchmark ---
    # return angle_matrix*composition[0] + angle_inclusion*composition[1]
    return angle_matrix

def sigma_yield(p_k, plastic_strain, composition):
    sigma_min = 1.0
    sigma_max = C0
    # value = p_k*sin(angle_friction(composition)) + cohesion(plastic_strain)*cos(angle_friction(composition))
    value = p_k*sin(friction_angle(plastic_strain)) + cohesion(plastic_strain)*cos(friction_angle(plastic_strain))

    return conditional(lt(value, sigma_min), sigma_min, value)
    # return conditional(lt(value, sigma_min), sigma_min, conditional(gt(value, sigma_max), sigma_max, value))

def plastic_strain_integration(mesh, tracers_in_cells, tracers, dt, yield_function, strain_rate, Temp):
    for j in range(mesh.num_cells()):
        for i in range(0,len(tracers_in_cells[j])):
            tracer_no = tracers_in_cells[j][i]
            px = tracers[tracer_no][0]
            py = tracers[tracer_no][1]
            temp_tracer = Temp(Point(px, py))

            yield_function_tracer = yield_function(Point(px,py))
            strain_rate_tracer = strain_rate(Point(px,py))

            if (yield_function_tracer >= -1e-3): # Add some toleration
                tracers[tracer_no][7] += float(dt)*strain_rate_tracer

            if (healing==True):
                tracers[tracer_no][7] = tracers[tracer_no][7]/(1.0+float(dt)/recovery_time)

# --- Elasticity ---

def z(visc, G, dt):
    return dt/(dt + visc/G)

def get_new_stress(mesh,Temp,xm,stress_invariant,strain_rate_invariant, p_k):
    ranks = []
    eval_type = "local"
    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()
        temp_cell = Temp(Point(centroid.x(),centroid.y()))
        eps_cell = strain_rate_invariant(Point(centroid.x(),centroid.y()))
        tau_cell = stress_invariant(Point(centroid.x(),centroid.y()))
        p_cell = p_k(Point(centroid.x(),centroid.y()))
        xm_cell = xm(Point(centroid.x(),centroid.y()))

        visc_cell = 1.0/(1.0/eta_diff(temp_cell)\
                        + 1.0/eta_disl(temp_cell,eps_cell,tau_cell,p_cell, eval_type) + 1.0/eta_max\
                        + 1.0/(eta_BS(temp_cell,eps_cell,tau_cell,p_cell)+eta_GBS(temp_cell,eps_cell,tau_cell,p_cell, eval_type)))*exp(-45.0*xm_cell)

        new_stress = 2.0*visc_cell*eps_cell + (1.0-z_cell)*tau_old_cell
        ranks.append(new_stress)
        
    ranks = numpy.array(ranks)                                                            
    stress_invariant.vector().set_local(ranks)

def get_new_stress_iter(mesh, Temp, xm, dt, stress_invariant, old_stress_invariant, strain_rate_invariant, step, iter, p_k):
    ranks = []
    eval_type = "local"

    for j in range(mesh.num_cells()):
        centroid = Cell(mesh, j).midpoint()
        temp_cell = Temp(Point(centroid.x(),centroid.y()))
        eps_cell = strain_rate_invariant(Point(centroid.x(),centroid.y()))
        tau_cell = stress_invariant(Point(centroid.x(),centroid.y()))
        tau_old_cell = old_stress_invariant(Point(centroid.x(),centroid.y()))
        p_cell = p_k(Point(centroid.x(),centroid.y()))
        xm_cell = xm(Point(centroid.x(),centroid.y()))
        gs_cell = d_grain(Point(centroid.x(),centroid.y()))
        G_cell = 3.52e9

        # --- First estimate of viscosity from previous stress ---
        if (step==0 or (step==1 and iter ==0) and (reloading_HDF5==False or reloading_tracers==False)):
            visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
        else:
            visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)\
                            + 1.0/eta_disl(temp_cell,eps_cell,tau_cell,p_cell, eval_type, gs_cell) + 1.0/eta_max\
                            + 1.0/(eta_BS(temp_cell,eps_cell,tau_cell,p_cell, gs_cell)+eta_GBS(temp_cell,eps_cell,tau_cell,p_cell, eval_type, gs_cell)))*exp(-45.0*xm_cell)

        stress_l = tau_old_cell
        stress_r = tau_old_cell

        if (elasticity==False):
            new_stress = 2.0*visc_cell*eps_cell
        else:
            z_cell = float(dt)/(float(dt) + visc_cell/G_cell)
            new_stress = 2.0*z_cell*visc_cell*eps_cell+ (1.0-z_cell)*tau_old_cell

        if (new_stress < tau_cell):
            stress_l = new_stress
        if (new_stress > tau_cell):
            stress_r = new_stress

        # --- Second estimate of viscosity from the new stress ---
        if (step==0 or (step==1 and iter ==0)):
            visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
        else:
            visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)\
                        + 1.0/eta_disl(temp_cell,eps_cell,new_stress,p_cell, eval_type, gs_cell) + 1.0/eta_max\
                        + 1.0/(eta_BS(temp_cell,eps_cell,new_stress,p_cell, gs_cell)+eta_GBS(temp_cell,eps_cell,new_stress,p_cell, eval_type, gs_cell)))*exp(-45.0*xm_cell)

        if (elasticity==False):
            new_stress = 2.0*visc_cell*eps_cell
        else:
            z_cell = float(dt)/(float(dt) + visc_cell/G_cell)
            new_stress = 2.0*z_cell*visc_cell*eps_cell+ (1.0-z_cell)*tau_old_cell

        if (new_stress < stress_l):
            stress_l = new_stress
        if (new_stress > stress_r):
            stress_r = new_stress

        # --- We have left and right boundary, let's iterate ---
        error = 1
        k=0
        # for k in range(0,iter_max):
        while (error > stress_iter_error and k < 20):
            stress_c = (stress_l+stress_r)/2.0

            if (step==0 or (step==1 and iter==0)):
                visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)  +  1.0/eta_max)*exp(-45.0*xm_cell)
            else:
                visc_cell = 1.0/(1.0/eta_diff(temp_cell, gs_cell)\
                            + 1.0/eta_disl(temp_cell,eps_cell,stress_c,p_cell, eval_type, gs_cell) + 1.0/eta_max\
                            + 1.0/(eta_BS(temp_cell,eps_cell,stress_c,p_cell, gs_cell)+eta_GBS(temp_cell,eps_cell,stress_c,p_cell, eval_type, gs_cell)))*exp(-45.0*xm_cell)

            if (elasticity==False):
                new_stress = 2.0*visc_cell*eps_cell
            else:
                z_cell = float(dt)/(float(dt) + visc_cell/G_cell)
                new_stress = 2.0*z_cell*visc_cell*eps_cell+ (1.0-z_cell)*tau_old_cell

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
    stress_invariant.vector().set_local(ranks)

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
        visc_disl = eta_disl(temp_cell, eps_cell, tau_cell, p_cell, eval_type)
        visc_GBS_BS = eta_GBS(temp_cell, eps_cell, tau_cell, p_cell, eval_type) + eta_BS(temp_cell, eps_cell, tau_cell, p_cell)
        
        if (plasticity==True):
            visc_plast = 0.5*ys_cell/max_function(eps_cell - (ys_cell - tau_old_cell)/(2.0*shear_modulus*float(dt)), float(sr_min)*1e-6)
            array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, visc_plast],[5.0, eta_max]]
        else:
            array_mech = [[1.0, visc_diff],[2.0, visc_GBS_BS],[3.0, visc_disl],[4.0, eta_max]]

        array_mech.sort(key=lambda x: x[1])

        ranks.append(array_mech[0][0])
    
    ranks = numpy.array(ranks)
    mechanisms.vector().set_local(ranks)


def eta_eff(p_k, strain_rate, Temp, xm, plastic_strain, stress, old_stress_invariant, step, dt, sr_min, composition, z_function):
    eval_type = "mesh"
    # --- Plastic viscosity ---
    if (plasticity == True):
        if ((elasticity==False or (elasticity==True and step==0)) and (reloading_HDF5==False or reloading_tracers==False)):
            eta_p = 0.5*sigma_yield(p_k,plastic_strain, composition)/strain_rate
        else:
            # eta_p = 0.5*sigma_yield(p_k,plastic_strain, composition)/max_function(strain_rate\
            #         - (sigma_yield(p_k,plastic_strain, composition) - old_stress_invariant)/(2.0*G*dt), sr_min*1e-6)
            
            # eta_p = 0.5*sigma_yield(p_k,plastic_strain, composition)/max_function(strain_rate\
            #         - (sigma_yield(p_k,plastic_strain, composition) - old_stress_invariant)/(2.0*G*dt), sr_min*1e-6)        
            
            eta_p = (sigma_yield(p_k,plastic_strain, composition) - (1-z_function)*old_stress_invariant)/(2.0*strain_rate*z_function)

    # --- Ductile viscosity ---
    if (viscosity_type=="constant"):
        eta_v = eta_0

    if (viscosity_type=="temp-dep"):
        # eta_v = eta_0*exp(60e3/R*(1.0/Temp - 1.0/T_bot))#*exp(-45.0*xm)
        # eta_v = eta_0*exp(-(Q_activ*(T_bot - T_top)/8.314/T_melt**2)*(Temp - T_bot)/170.0)*exp(-45.0*xm)
        # eta_v = eta_0*exp(-(Q_activ*(Temp - T_melt)/8.314/T_melt**2))*exp(-45.0*xm)
        # eta_v = eta_0/exp(Q_activ/8.314/T_melt)*exp(Q_activ/8.314/Temp)*exp(-45.0*xm)
        eta_v = eta_0*exp(-14.0*(Temp - T_bot)/170.0)*exp(-45.0*xm)

    if (viscosity_type=="GK_2001"):
        if ((step == 0 or (step == 1 and iter == 0)) and (reloading_HDF5 == False or reloading_tracers == False)):
            eta_v = 1.0/(1.0/eta_diff(Temp) + 1.0/eta_max)*exp(-45.0*xm) 

        else:
            eta_v = 1.0/(1.0/eta_diff(Temp)\
                    + 1.0/eta_disl(Temp ,strain_rate, stress, p_k, eval_type)\
                    + 1.0/(eta_GBS(Temp, strain_rate, stress, p_k, eval_type) + eta_BS(Temp, strain_rate, stress, p_k))\
                    + 1.0/eta_max)*exp(-45.0*xm)
                        
    if (viscosity_type == "composition"):
        # --- Dense layer benchmark ---
        # eta_v = 1.0e0**composition[0] * 1.0**composition[1]\

        # --- Rising plume benchmark ---
        # eta_v = eta_mantle**composition[0] * eta_lid**composition[1] * eta_plume**composition[2]  

        # --- Shear bands benchmark ---
        eta_v = 1e25**composition[0] * 1.0e20**composition[1]

    # --- NOT AVERAGING the viscosities ---
    if (plasticity==True):
        if (step==0 or (step==1 and iter ==0)):
            return eta_v
        else:
            return conditional(lt(eta_p, eta_min_plast), eta_min_plast, conditional(lt(eta_p, eta_v), eta_p, eta_v))
    else:
        return conditional(lt(eta_v, eta_max), eta_v, eta_max)

def shear_modulus(composition):
    # --- Rising plume benchmark ---
    # return G_mantle**composition[0] * G_lid**composition[1] * G_plume**composition[2]

    return 5e10

def eta(Q, A, n, m, Temp, strain_rate_II, stress_II, pressure):
    # See Gerya (2009) chapter 6.2 for derivation of prefactors
    # The correction applies also to diffusion creep, see the function below and eq. (6.8b) in Gerya (2009)
    if (viscosity_from=="strain_rate"):
        return 1.0/(3.0**((n+1)/(2.0*n))*2.0**((n-1.0)/n))*A**(-1.0/n)*d_grain**(m/n)*strain_rate_II**((1.0-n)/n)*exp((Q+pressure*V_act)/(n*R*Temp))
        
    if (viscosity_from=="stress"):
        return 1.0/(3.0**((n+1)/2.0))*A**(-1.0)*d_grain**m*stress_II**(1.0-n)*exp((Q+pressure*V_act)/(R*Temp))

def eta_diff(Temp):
    pi=3.1415926535
    V_m = 1.97e-5           #m^3
    D_0v = 9.10e-4          #m^2 s^-1
    D_0b = 6.4e-4           #m^2 s^-1
    Q_v = 59.4e3            #J mol^-1
    Q_b = 49.0e3            #J mol^-1
    delta= 9.04e-10         #m

    D_v = D_0v*exp(-Q_v/(R*Temp))
    D_b = D_0b*exp(-Q_b/(R*Temp))

    # 3/2 come from the scaling between diff. stress, see Eq. (6.8b) in Gerya (2009)
    return 1.0/(3.0/2.0*84*V_m)*R*Temp*d_grain**2/(D_v + pi*delta*D_b/d_grain)

def eta_disl(Temp, strain_rate_II, stress_II, pressure, eval_type):
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
        eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II, pressure),\
        eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II, pressure))

    if (eval_type == "local"):
        if (Temp > T_crit):
            return eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II, pressure)
        else:
            return eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II, pressure)

def eta_GBS(Temp, strain_rate_II, stress_II, pressure, eval_type):
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
        eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II, pressure),\
        eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II, pressure))  

    if (eval_type == "local"):
        if (Temp > T_crit):
            return eta(Q_above, A_above, n, p, Temp, strain_rate_II, stress_II, pressure)
        else:
            return eta(Q_below, A_below, n, p, Temp, strain_rate_II, stress_II, pressure)    

def eta_BS(Temp,strain_rate_II,stress_II,pressure):
        Q = 60e3          #J/mol
        A = 2.2e-7        #Pa^{-2.4} s^{-1}
        n = 2.4           # -
        p = 0.0           # -

        return eta(Q,A,n,p,Temp,strain_rate_II,stress_II,pressure)
