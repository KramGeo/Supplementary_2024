# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 17:06:30 2021

@author: ocollet001
"""
import numpy as np

#%% HELPER FUNCTIONS
""" FUNCTION poissonRatioFromModuli """
def poissonRatioFromModuli(K,G):
    """
    Calculate Poisson's ratio from input bulk and shear moduli.
    
    Parameters
    ----------
    K : float
        bulk modulus
    G : float
        shear modulus
    
    Notes
    -----
    Both inputs must have the same units.
    """
    poisson = (1.5*K-G)/(3.0*K+G);
    return poisson

#%% GRANULAR MODELS

""" FUNCTION HertzMindlin """
def HertzMindlin(KmGPa, GmGPa, PeffMPa,
                 phi=0.40, C=9., RR=1., FC=1.):
    """
    Calculate bulk and shear moduli using Hertz-Mindlin theory.

    Parameters
    ----------
    KmGPa : float
        Bulk modulus of solid (grain) (in GPa).
    GmGPa : float
        Shear modulus of solid (grain) (in GPa).
    PeffMPa : float
        Effective pressure (in MPa).
    phi : float, optional
        Porosity. The default is 0.40.
    C : float, optional
        Coordination number. The default is 9.
    RR : float, optional
        Grain angularity parameter. The default is 1.
    FC : float, optional
        Friction coefficient. The default is 1.

    Returns
    -------
    Khm : float
        Effective bulk modulus (in GPa).
    Ghm : float
        Effective shear modulus (in GPa).

    """             
    PRm = poissonRatioFromModuli(KmGPa,GmGPa);   
    
    PeffGpa = PeffMPa/1000.;   
    
    chi = C*(1.-phi)*GmGPa/(np.pi*(1.-PRm));
    chi2 = chi*chi;    
    
    Khm = ( (RR*chi2*PeffGpa) / 18.0 ) ** (1./3.);
    #Khm = ( (C*C*RR*(1.-Phic)*(1.-Phic)*GmGPa*GmGPa*(PeffMPa/1000.)) / (18*9.8696*(1-PRm)*(1-PRm)) ) ** (1./3.);

    Ghm = (2.+3.*FC-PRm*(1.+3.*FC)) / (10.-5.*PRm) * ( 1.5*RR*chi2*PeffGpa ) ** (1./3.);
    #Ghm = (2.+3.*FC-PRm*(1.+3.*FC)) / (10.-5.*PRm) * ( (3.*C*C*RR*(1-Phic)*(1-Phic)*GmGPa*GmGPa*PeffGpa) / (2.*9.8696*(1.-PRm)*(1-PRm)) ) ** (1./3.);

    return Khm, Ghm

""" FUNCTION Walton """
def Walton(KmGPa, GmGPa, PeffMPa, phi=0.40, C=9, R=3/5):
    """
    Calculate bulk and shear moduli using Walton (1987) model.
    See eq. (9.49) Pride, Hydrogeophysics (2005).

    Parameters
    ----------
    KmGPa : float
        Bulk modulus of solid (grain) (in GPa).
    GmGPa : float
        Shear modulus of solid (grain) (in GPa).
    PeffMPa : float
        Effective pressure (in MPa).
    phi : float, optional
        Porosity. The default is 0.40.
    C : float, optional
        Coordination number. The default is 9.
    R : float, optional
        Shear to bulk modulus ratio. The default is 3/5.

    Returns
    -------
    Kd : float
        Effective drained bulk modulus (in GPa).
    Gd : float
        Effective drained shear modulus (in GPa).

    """   
    PeffGPa = PeffMPa/1000.;
    Cm = (1/(4*np.pi))*(1/GmGPa + 1/(KmGPa+GmGPa/3))
           
    Kd = (1/6)*(3*(1-phi)**2*C**2*PeffGPa/(np.pi**4*Cm**2))**(1/3)
    Gd = R*Kd

    return Kd, Gd

""" FUNCTION Walton_modified """
def Walton_modified(KmGPa, GmGPa, PeffMPa, P0MPa, phi=0.40, C0=9, R=3/5):
    """
    Calculate bulk and shear moduli using a modified version of Walton (1987)
    model, which accounts for the negative dilatation occuring at low pressures.
    See eq. (9.53) Pride, Hydrogeophysics (2005).

    Parameters
    ----------
    KmGPa : float
        Bulk modulus of solid (grain) (in GPa).
    GmGPa : float
        Shear modulus of solid (grain) (in GPa).
    PeffMPa : float
        Effective pressure (in MPa).
    P0MPa : float
        Critical pressure for which C=C0 (in MPa).
    phi : float, optional
        Porosity. The default is 0.40.
    C0 : float, optional
        Coordination number in the high-stress limit. The default is 9.
    R : float, optional
        Shear to bulk modulus ratio. The default is 3/5.

    Returns
    -------
    Kd : float
        Effective drained bulk modulus (in GPa).
    Gd : float
        Effective drained shear modulus (in GPa).

    """   
    PeffGPa = PeffMPa/1000.;
    P0GPa = P0MPa/1000.;
    Cm = (1/(4*np.pi))*(1/GmGPa + 1/(KmGPa+GmGPa/3))
           
    Kd = (1/6)*(4*(1-phi)**2*C0**2*P0GPa/(np.pi**4*Cm**2))**(1/3)*(PeffGPa/P0GPa)**(1/2)/(1+(16*PeffGPa/(9*P0GPa))**4)**(1/24)
    Gd = R*Kd

    return Kd, Gd