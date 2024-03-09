# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 10:45:32 2022

@author: ocollet001
"""
import numpy as np

#%% Van Genuchten model

def Swe_from_wt_VG(z, wt, alpha_vg, n_vg):
    """
    Compute effective water saturation using Van Genuchten (VG) model.
    See eq. 16 Solazzi et al. (2021) 

    Parameters
    ----------
    z : numpy array
        Depth (m).
    wt : float
        Water table depth (m).
    alpha_vg : float
        Inverse of the entry pressure (1/m).
    n_vg : float
        Pore size distribution parameter (unitless).

    Returns
    -------
    Effective water saturation (float or numpy array).

    """
    m_vg = 1-1/n_vg
    
    # Initialize output
    Swe = np.ones(np.shape(z))
    
    # Take absolute values if depths are counted negative from surface
    z = np.abs(z)
    wt = np.abs(wt)
    
    # Compute effective saturation above water table using VG model
    I = np.where(z<wt)
    Swe[I] = (1+(alpha_vg*(wt-z[I]))**n_vg)**(-m_vg)
    
    return Swe
    