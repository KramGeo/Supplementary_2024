# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 11:24:56 2021

@author: ocollet001
"""
import numpy as np

def velocities_polar_angle(C11,C33,C13,C44,C66,rho,theta):
    """
    Calculate velocities as a function of angle in TI anisotropic media.
    The symmetry axis is assumed to be along x3.
    See e.g. Thomsen - SEG DISC 2002 - p.1-28.   

    Parameters
    ----------
    C11 : float
        C11 stiffness component (GPa).
    C33 : float
        C33 stiffness component (GPa).
    C13 : float
        C13 stiffness component (GPa).
    C44 : float
        C44 stiffness component (GPa).
    C66 : float
        C66 stiffness component (GPa).
    rho : float
        Density (g/cc).
    theta : float array
        Angle between wave propagation and symmetry axis (rad).

    Returns
    -------
    Vp, Vs_perp, Vs_par : float arrays
        P_wave, S_perp (i.e. SV) and S_par (i.e. SH) wave velocities (km/s).

    """
    cos2theta = np.cos(theta)*np.cos(theta)
    sin2theta = np.sin(theta)*np.sin(theta)
    sin4theta = sin2theta*sin2theta
    
    D = np.sqrt((C33-C44)**2 +
         2*(2*(C13+C44)**2-(C33-C44)*(C11+C33-2*C44))*sin2theta +
         ((C11+C33-2*C44)**2-4*(C13+C44)**2)*sin4theta)
    
    Vp = np.sqrt((C33+C44+(C11-C33)*sin2theta+D)/(2*rho))
    Vs_perp = np.sqrt((C33+C44+(C11-C33)*sin2theta-D)/(2*rho))
    Vs_par = np.sqrt((C44*cos2theta+C66*sin2theta)/rho)

    return Vp, Vs_perp, Vs_par
    