# -*- coding: utf-8 -*-
"""
bounds2.py 
calculate bounds using Voight, Reuss, Hill and Hashin-Shritkman for 2 phase mineral
   Hashin-Shritkman based on hashin_shtrikman_bounds.m
   see rock physics handbook, follows Hasin-Shtirkman-Wapole form for 2 phases
   or see Dvorkin's book eqn 2.7  (most straight forward to follow)
Created on Tue Jan 09 11:15:32 2018
@author: jdownton
"""

def VoightAvr(M1,f1,M2):
    M=f1*M1+(1-f1)*M2
    return M    
    
def ReussAvr(M1,f1,M2):
    M=1/(f1/M1+(1-f1)/M2)
    return M

def HillAvr(M1,f1,M2, weight=0.5):
    M = weight*VoightAvr(M1,f1,M2) + (1-weight)*ReussAvr(M1,f1,M2)
    return M

#calculate the Hashin Shtrikman upper (/lower) bound for bulk moduli, this is a helper routine for hashin_shtrikman_bounds()
def bulk_HS_bound(z, K1, K2, f1, f2):
    z43 = 4*z/3;
    return 1./(f1/(K1 + z43) + f2/(K2 + z43)) - z43;

#calculate the Hashin Shtrikman upper (/lower) bound for shear moduli, this is a helper routine for hashin_shtrikman_bounds()
def shear_HS_bound(z, G1, G2, f1, f2):
    denom=(f1*(G2 + z) + f2*(G1 + z));
    if (denom <= 0):
        gam=G1;
    else:
        gam=(G1 + z)*(G2 + z)/denom - z;
    return gam;

#calculate z used in the calculation of the Hashin Shtrikman upper (/lower) bound for shear moduli, this is a helper routine for hashin_shtrikman_bounds()
def z_HS(K, G):
    return G/6.*(9.*K + 8.*G)/(K + 2.*G);  #since K>0 and G>=0 the denominator should never be zero

#Returns Hashin-Shtrikman upper and lower bounds on bulk and shear modulus.
# 
def hashin_shtrikman(K1,  #Bulk modulus of material 1 <units=GPa>
                    G1,  #Shear modulus of material 1 <units=GPa>
                    f1,  #Material 1 fraction <units=V/V>
                    K2,  #Bulk modulus of material 2 <units=GPa>
                    G2): #Shear modulus of material 2 <units=GPa>

  f1 = max(0, min(f1, 1));  #Material 1 fraction <units=V/V>
  f2 = 1 - f1;              #Material 2 fraction <units=V/V>

  Gmax = max(G1, G2);
  Gmin = min(G1, G2);
  Kmax = max(K1, K2);
  Kmin = min(K1, K2);

  Khsup = bulk_HS_bound(Gmax, K1, K2, f1, f2);
  Khslo = bulk_HS_bound(Gmin, K1, K2, f1, f2);

  Zmax = z_HS(Kmax, Gmax);
  Zmin = z_HS(Kmin, Gmin);

  Ghsup = shear_HS_bound(Zmax, G1, G2, f1, f2);
  Ghslo = shear_HS_bound(Zmin, G1, G2, f1, f2);

  return Khsup, Ghsup, Khslo, Ghslo;