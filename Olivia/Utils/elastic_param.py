# -*- coding: utf-8 -*-
"""
Utilities to compute elastic parameters.

Isotropic media:
    K : bulk modulus
    G : shear modulus
    E : Young's modulus
    nu: Poissons'ratio
    Vp: P-wave velocity
    Vs: S-wave velocity

Created on Fri Jun 25 13:32:59 2021

@author: ocollet001
"""
import numpy as np

def K_from_E_nu(E, nu):
    """
    Return bulk modulus from Young's modulus and Poisson's ratio.

    Parameters
    ----------
    E : float array
        Young's modulus (GPa).
    nu : float array
        Poisson's ratio.

    Returns
    -------
    K : float array
        Bulk modulus (GPa).

    """
    return E/(3*(1-2*nu))

def G_from_E_nu(E, nu):
    """
    Return shear modulus from Young's modulus and Poisson's ratio.

    Parameters
    ----------
    E : float array
        Young's modulus (GPa).
    nu : float array
        Poisson's ratio.

    Returns
    -------
    G : float array
        Shear modulus (GPa).

    """
    return E/(2*(1+nu))

def nu_from_K_G(K, G):
    """
    Return Poisson's ratio from bulk and shear moduli.

    Parameters
    ----------
    K : float array
        Bulk modulus (GPa).
    G : float array
        Shear modulus (GPa).

    Returns
    -------
    nu : float array
        Poisson's ratio (unitless).

    """
    return 0.5*(3*K-2*G)/(3*K+G)

def Vp_from_K_G(K,G,rho):
    """
    Calculate P-wave velocity from bulk and shear moduli and density.

    Parameters
    ----------
    K : float array
        Bulk modulus (GPa).
    G : float array
        Shear modulus (GPa).
    rho : float array
        Density (g/cc).    

    Returns
    -------
    Vp : float array
        P-wave velocity (km/s).

    """
    return np.sqrt((K+4*G/3)/rho)

def Vs_from_G(G,rho):
    """
    Calculate S-wave velocity from shear modulus and density.

    Parameters
    ----------
    G : float array
        Shear modulus (GPa).
    rho : float array
        Density (g/cc).    

    Returns
    -------
    Vs : float array
        S-wave velocity (km/s).

    """
    return np.sqrt(G/rho)

def TI_Young_Poisson_from_stiffness(C11,C33,C13,C44,C66):
    """
    Young's moduli and Poisson's ratio for TI media
    with symmetry axis along x3.
    See Mavko's book p.36

    Parameters
    ----------
    C11 : TYPE
        DESCRIPTION.
    C33 : TYPE
        DESCRIPTION.
    C13 : TYPE
        DESCRIPTION.
    C44 : TYPE
        DESCRIPTION.
    C66 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    C12 = C11-2*C66
    
    E33 = C33 - (2*C13*C13)/(C11+C12)
    nu31 = C13/(C11+C12)
    
    E11 = C11 + (C13*C13*(C12-C11)+C12*(C13*C13-C33*C12))/(C33*C11-C13*C13)
    nu13 = (C13*(C11-C12))/(C33*C11-C13*C13)
    nu12 = (C33*C12-C13*C13)/(C33*C11-C13*C13)
    
    return E11, E33, nu13, nu12, nu31