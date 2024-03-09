# -*- coding: utf-8 -*-
"""
Fucntions to compute Voigt, Reuss and Hill averages.
Created on Tue Dec  7 09:18:28 2021

@author: ocollet001
"""

def VoigtAvr(M1,f1,M2):
    """
    Calculate Voigt average (two phase).

    Parameters
    ----------
    M1 : float or float array
        Elastic modulus of material 1.
    f1 : float or float array
        Volume fraction of material 1.
    M2 : float or float array
        Elastic modulus of material 2.

    Returns
    -------
    M : float or float array
        Effective elastic modulus computed using Voigt average.

    """
    M = f1*M1+(1-f1)*M2
    return M    
    
def ReussAvr(M1,f1,M2):
    """
    Calculate Reuss average (two phase).

    Parameters
    ----------
    M1 : float or float array
        Elastic modulus of material 1.
    f1 : float or float array
        Volume fraction of material 1.
    M2 : float or float array
        Elastic modulus of material 2.

    Returns
    -------
    M : float or float array
        Effective elastic modulus computed using Reuss average.

    """
    M = 1/(f1/M1+(1-f1)/M2)
    return M

def HillAvr(M1,f1,M2, weight=0.5):
    """
    Calculate Hill average (two phase).

    Parameters
    ----------
    M1 : float or float array
        Elastic modulus of material 1.
    f1 : float or float array
        Volume fraction of material 1.
    M2 : float or float array
        Elastic modulus of material 2.
    weight : float, default=0.5
        Weight to be applied to Voigt average. 
        If weight=0.5: standard Hill average;
        Else: modified Hill's average.

    Returns
    -------
    M : float or float array
        Effective elastic modulus computed using Reuss average.

    """
    M = weight*VoigtAvr(M1,f1,M2) + (1-weight)*ReussAvr(M1,f1,M2)
    return M