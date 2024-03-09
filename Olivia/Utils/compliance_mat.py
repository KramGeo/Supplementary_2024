# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:26:30 2021

@author: ocollet001
"""
import numpy as np

def iso_compliance_mat(E, nu):
    """
    Return compliance matrix of isotropic media

    Parameters
    ----------
    E : float
        Young's modulus.
    nu : float
        Poisson's ratio.

    Returns
    -------
    Compliance matrix (6x6 numpy array).
    
    """
    test = False
    if test:
        mu = 2*E/(1+nu)
    else:
        mu = E/(2*(1+nu))
    
    S0 = np.array([[1/E,   -nu/E, -nu/E,  0.,    0.,    0.],
                   [-nu/E,  1/E,  -nu/E,  0.,    0.,    0.],
                   [-nu/E, -nu/E,  1/E,   0.,    0.,    0.],
                   [0.,    0.,    0.,   1/mu,    0.,    0.],
                   [0.,    0.,    0.,     0.,    1/mu,  0.],
                   [0.,    0.,    0.,     0.,    0.,    1/mu]])
    
    return S0
