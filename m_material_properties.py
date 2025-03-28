from dolfin import *
import numpy
import os 
from m_parameters import *

def density(Temp, composition, xm):
        """Evaluates the density :math:`\\rho`\ . 

        :param Temp: temperature (:math:`T`\ )
        :param composition: material composition (:math:`C`\ )
        :param xm: melt fraction (:math:`\chi`\ )

        :returns:
           * For ice :math:`\\rho = (1 - \\chi)\\rho_i(T) + \\chi\\rho_w`
           * Another user-defined function of temperature, composition and melt fraction. 
             Formulas used in the benchmarks are predefined and commented.

        """
        # Density of ice following Rottger et al. (1994) and Feistel and Wagner (2006)
        a0 = 128.2147
        a3 = -1.3152e-6
        a4 = 2.4837e-8
        a5 = -1.6064e-10
        a6 = 4.6097e-13
        a7 = -4.966e-16

        VV = a0 + a3*Temp**3 + a4*Temp**4  + a5*Temp**5 + a6*Temp**6 + a7*Temp**7
        mm = a0 + a3*T_ref**3 + a4*T_ref**4  + a5*T_ref**5 + a6*T_ref**6 + a7*T_ref**7

        # --- Van Keken RT-instability benchmark ---
        # return -1*composition[0] 

        # --- Shear bands benchmark ---
        # return rho_s

        # --- Ice with melt ---
        return (1.0-xm)*rho_s*mm/VV + xm*rho_m

        # --- Ice ---
        if (len(materials) == 0):
                # With melt
                return (1.0-xm)*rho_s*mm/VV + xm*rho_m
        else:
                # With salt and melt
                return ((1.0-xm)*rho_s*mm/VV + xm*rho_m)*(composition[0] + 1.005*composition[1])

        # return rho_s*(1.0-alpha_exp*(Temp-T_ref)) + xm*rho_s*(1.0-rho_s/rho_l) # Tobie et al. (2003)
        # return (1.0-(xm+5e-2))*rho_s*(1.0-alpha_exp*(Temp-T_ref)) + (xm+5e-2)*rho_m # Klara PhD. (2015)
        #     return -Ra*Temp - Rb*composition[0] 

    # return Ra/(2.5e-5*1e3) - Ra*Temp
#     return rho_mantle*composition[0] + rho_lid*composition[1] + rho_plume*composition[2]
        # return 2700.0
        
def k(Temp, composition):
        """Evaluates the thermal conductivity :math:`k`\ . 

        :param Temp: temperature (:math:`T`\ )
        :param composition: material composition (:math:`C`\ )

        :returns: A user-defined function of temperature and composition.

        """
        return 2.3 #567.0/Temp

def cp(Temp, composition):
        """Evaluates the specific heat capacity at constant pressure :math:`c_p`\ . 

        :param Temp: temperature (:math:`T`\ )
        :param composition: material composition (:math:`C`\ )

        :returns: A user-defined function of temperature and composition.

        """
        return 185.0 + 7.037*Temp

# def shear_modulus(composition):
#     # --- Rising plume benchmark ---
#     # return G_mantle**composition[0] * G_lid**composition[1] * G_plume**composition[2]

#     return 5e10