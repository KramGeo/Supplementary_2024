# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 13:49:45 2021

@author: ocollet001
"""

def ani_param_HTI(C11, C33, C13, C44, C55, 
                  delta_ap=False, eta_ap=False):
    """
    Return anisotropy parameters for HTI media 
    (symmetry axis is along x1).

    Parameters
    ----------
    C11 : float array
        C11 stiffness component.
    C33 : float array
        C33 stiffness component.
    C13 : float array
        C13 stiffness component.
    C44 : float array
        C44 stiffness component.
    C55 : float array
        C55 stiffness component.
    delta_ap : bool, default=False
        Flag to control whether weak anisotropy approx is used for delta.    
    eta_ap : bool, default=False
        Flag to control whether weak anisotropy approx is used for eta.
        
    Returns
    -------
    eps, delta, gamma, eta : float arrays
        Anisotropy parameters for HTI media.

    """
    eps = 0.5*(C33-C11)/C11
    
    if (delta_ap): #weak anisotropy approximation
        delta = (C13+2*C55-C11)/C11
    else:  #exact expression (default)
        delta = 0.5*((C13+C55)**2-(C11-C55)**2)/(C11*(C11-C55))
        
    gamma = 0.5*(C44-C55)/C55
    
    if (eta_ap): #weak anisotropy approximation
        eta = eps-delta
    else: #exact expression (default)
        eta = (eps-delta)/(1+2*delta)
        
    return eps, delta, gamma, eta

def ani_param_VTI(C11, C33, C13, C44, C66, 
                  delta_ap=False, eta_ap=False):
    """
    Return anisotropy parameters for HTI media 
    (symmetry axis is along x3).

    Parameters
    ----------
    C11 : float array
        C11 stiffness component.
    C33 : float array
        C33 stiffness component.
    C13 : float array
        C13 stiffness component.
    C44 : float array
        C44 stiffness component.
    C66 : float array
        C66 stiffness component.
    delta_ap : bool, default=False
        Flag to control whether weak anisotropy approx is used for delta.    
    eta_ap : bool, default=False
        Flag to control whether weak anisotropy approx is used for eta.
        
    Returns
    -------
    eps, delta, gamma, eta : float arrays
        Anisotropy parameters for VTI media.

    """
    eps = 0.5*(C11-C33)/C33
    
    if (delta_ap): #weak anisotropy approximation
        delta = (C13+2*C44-C33)/C33
    else:  #exact expression (default)
        delta = 0.5*((C13+C44)**2-(C33-C44)**2)/(C33*(C33-C44))
        
    gamma = 0.5*(C66-C44)/C44
    
    if (eta_ap): #weak anisotropy approximation
        eta = eps-delta
    else: #exact expression (default)
        eta = (eps-delta)/(1+2*delta)
        
    return eps, delta, gamma, eta