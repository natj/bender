import numpy as np


constbb = 5.040366e22  # hzkeV^4 * 2*h/c^2
ergkev  = 6.24150934326e8 # erg/keV


##################################################
# Planck function 

#Planck function B_E(T) = erg cm^-2 s^-1 str^-1 keV^-1
def BE(T, E):
    return constbb * (E**3)/np.expm1(E/T)

#Planck function photon flux N_E(T) = ph cm^-2 s^-1 str^-1 keV^-1
def NBE(T, E):
    return constbb * (E**3)/np.expm1(E/T)/E*ergkev

#Number of black body photons N(T) = \int dE B_E(T)/E(keV) * ergkev
def NB(T):
    return constbb*ergkev*2.404*T**3

#Energy of black body photons E(T) = \int dE B_E(T)
def EB(T):
    return constbb*ergkev*(np.pi**4/15.0)*T**4


##################################################
#Beaming, i.e. angle dependency of the radiation

#Isotropic beaming
def isotropic_beaming(mu):
    return 1.0






