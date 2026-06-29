# --- Python modules ---
from dolfin import *
import numpy
import os 

# --- MilleFEuiIle modules ---
from m_parameters import *

def rho_ice_water(Temp, xm):
        a0 = 128.2147
        a3 = -1.3152e-6
        a4 = 2.4837e-8
        a5 = -1.6064e-10
        a6 = 4.6097e-13
        a7 = -4.966e-16

        VV = a0 + a3*Temp**3 + a4*Temp**4  + a5*Temp**5 + a6*Temp**6 + a7*Temp**7
        mm = a0 + a3*T_ref**3 + a4*T_ref**4  + a5*T_ref**5 + a6*T_ref**6 + a7*T_ref**7

        return (1.0-xm)*rho_s*mm/VV + xm*rho_m


def rho(Temp, composition, xm):
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

        # --- Default ice with melt ---
        return (1.0-xm)*rho_s*mm/VV + xm*rho_m

        # --- Demo 4 ---
        # return ((1.0-xm)*rho_s*mm/VV + xm*rho_m)*(composition[0] + 1.004*composition[1])

def k(Temp, composition):
        """Evaluates the thermal conductivity :math:`k`\ . 

        :param Temp: temperature (:math:`T`\ )
        :param composition: material composition (:math:`C`\ )

        :returns: A user-defined function of temperature and composition.

        """
        k_H2O = 612.0/Temp
        
        # --- Demo 5 ---
        # k_H2O = 2.3

        return k_H2O

def cp(Temp, composition):
        """Evaluates the specific heat capacity at constant pressure :math:`c_p`\ . 

        :param Temp: temperature (:math:`T`\ )
        :param composition: material composition (:math:`C`\ )

        :returns: A user-defined function of temperature and composition.

        """
        cp_H2O = 185.0 + 7.037*Temp
        
        # --- Demo 5 ---
        # cp_H2O = 2100
        
        return cp_H2O