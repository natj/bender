

##################################################
#Physical conversion factors 

#Particle Data Group (Review of particle Physics 2015) values
# in SI units
c      = 2.99792458e8 #speed of light
h      = 6.62607004081e-34 #Planck constant
G      = 6.6730831e-11 #gravitational constant
eV     = 1.602176620898e-19 #electron volt/Joule
parsec = 3.08567758149e16 #parsec in m
Rs     = 2.9532500772e3 #Sch radius of Sun in m




##################################################
# Operate in units where G = c = 1.
# Use units of solar masses

##################################################
# (G Msun /c^2 )^-1
solar_mass_per_km = 2.0e3/Rs

# (G Msun/c^3)^-1
solar_mass_per_s = 2.0*c/Rs



##################################################
# other conversion factors
km_per_kpc = parsec #km / kpc = m/pc
cm_per_tenkpc = 1.0e-6 / parsec

#kelvin_per_kev = 1.16045e7

# Blackbody constant
# hzkeV^4 * 2*h/c^2
constbb = 1.0e15 * (eV/h)**4 * 2.0*h/c**2

erg_per_kev = 1.0e-10 / eV



# other constants
##################################################
#mass of sun
#Msun   = 1.988435e30
#kg_per_Msun    = 1.988435e30


