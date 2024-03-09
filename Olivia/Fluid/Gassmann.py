# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:58:48 2021

@author: ocollet001
"""

#%% FLUID MIXING LOWS

def fluidWood(S_f1, K_f1=2.25, K_f2=1e-4):
    """
    Calculate the effective fluid bulk modulus using Wood's equation.

    Parameters
    ----------
    S_f1 : float
        Saturation of the first fluid phase.
    K_f1 : float, optional
        Bulk modulus of the first fluid. The default is 2.25 GPa (water).
    K_f2 : TYPE, optional
        Bulk modulus of the second fluid. The default is 1e-4 GPa (air).

    Returns
    -------
    Kf : float
        Effective fluid bulk modulus (GPa).

    """
    
    Kf = 1/(S_f1/K_f1 + (1-S_f1)/K_f2)
    return Kf


#%% GASSMANN EQUATIONS

def isoGassmann(K,G,Kg,Kf,phi):
    """
    Calculate elastic param in saturated isotropic media 
    using the isotropic Gassmann equations.

    Parameters
    ----------
    K : float array
        Bulk modulus of the dry medium.
    G : float array
        Shear modulus of the dry medium.
    Kg : float
        Grain bulk modulus.
    Kf : float
        Fluid bulk modulus.
    phi : float
        Porosity.

    Returns
    -------
    Ksat : float array
        Bulk modulus of the saturated medium.
    
    Notes
    -----
    Units for K, G, Kg and Kf must be the same.
    Need logic if phi[i]=0.
    
    """
    alpha = 1-K/Kg
    M = Kg/((1-K/Kg)-phi*(1-Kg/Kf))
    Ksat = K + alpha*alpha*M
    
    return Ksat
    
def aniGassmann_HTI(C11,C33,C13,C44,C55,Kg,Kf,phi):
    """
    Calculate stiffness coeff in saturated HTI media 
    using the anisotropic Gassmann equations.
    No approximation is made.

    Parameters
    ----------
    C11, C33, C13, C44, C55 : float arrays
        Stiffness coefficients of the dry medium.
    Kg : float
        Grain bulk modulus.
    Kf : float
        Fluid bulk modulus.
    phi : float
        Porosity.

    Returns
    -------
    C11sat, C33sat, C13sat, C44sat, C55sat : float arrays
        Stiffness coefficients of the saturated medium.
    
    Notes
    -----
    Units for K and Cij must be the same.
    
    """
    C23 = C33-2*C44
    
    alpha1 = 1-(C11+2*C13)/(3*Kg)
    alpha2 = 1-(C13+C23+C33)/(3*Kg)
    alpha3 = alpha2
    alpha4 = 0
    alpha5 = 0

    Kstar = (1/9)*(C11+2*C33+4*C13+2*(C33-2*C44))
    M = Kg/((1-Kstar/Kg)-phi*(1-Kg/Kf))

    C11sat = C11 + alpha1*alpha1*M
    C13sat = C13 + alpha1*alpha3*M
    C33sat = C33 + alpha3*alpha3*M
    C44sat = C44 + alpha4*alpha4*M
    C55sat = C55 + alpha5*alpha5*M
    
    return C11sat, C33sat, C13sat, C44sat, C55sat