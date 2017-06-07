import numpy as np
import units


##################################################
# Planck function 

#Planck function B_E(T) = erg cm^-2 s^-1 str^-1 keV^-1
def BE(T, E):
    return units.constbb * (E**3)/np.expm1(E/T)

#Planck function photon flux N_E(T) = ph cm^-2 s^-1 str^-1 keV^-1
def NBE(T, E):
    return units.constbb * (E**3)/np.expm1(E/T)/E* units.erg_per_kev

#Number of black body photons N(T) = \int dE B_E(T)/E(keV) * ergkev
def NB(T):
    return units.constbb* units.erg_per_kev*2.404*T**3

#Energy of black body photons E(T) = \int dE B_E(T)
def EB(T):
    return units.constbb* units.erg_per_kev*(np.pi**4/15.0)*T**4


##################################################
#Beaming, i.e. angle dependency of the radiation

#Isotropic beaming
def isotropic_beaming(mu):
    return 1.0






