using Winston
using Colors
using Interpolations
using Cubature



######################
# Physical constants
const G = 6.67384e-8
const c = 2.99792458e10
const Msun = 1.9885469e33 #XXX
const km = 1.0e5
const ergkev = 6.24150934326e8 # erg/keV 
const cm_parsec = 3.24077929e-23 #1cm/10kpc # 3.08567758135
const constbb = 5.040366e22


#JP constants
#const G = 6.67384e-8 
#const c = 2.99792458e10  
#const Msun = 1.98892e33
#const km = 1.0e5
#const ergkev = 6.2415e8
#const cm_parsec =  3.2404e-23
#const constbb = 5.0396173e22 


#Physical configuration
incl = deg2rad(15.0)
M    = 1.6Msun
R    = 12.0km
fs   = 400


#ray tracing grid setup
Nrad = 100
Nchi = 100





#derived dimensionless values
const sini = sin(incl)
#const Rg = G*M/c^2
const Rg = 1.0
const X = G*M/(R*c^2)
const Osb = (2pi*fs)*sqrt(R^3/(G*M))
const U = 2*G*M/(R*c^2)
const imgscale = (G*M/c^2)^2

println("x=$X ($U) and Osb=$Osb incl=$incl")

#Hartle-Thorne parameters
const beta = 0.4454*Osb^2*X #Matter quadrupole moment; AlGendy & Morsink 2014
const quad = -0.11*(Osb/X)^2 #Energy quadrupole moment; AlGendy & Morsink 2014


I = sqrt(X)*(1.136 - 2.53*X + 5.6*X^2)*M*R^2 #Moment of inertia; AlGendy & Morsink 2014
J = I*(2pi*fs)


#dimensionless inertia
imom = sqrt(X)*(1.136 - 2.53*X + 5.6*X^2) #dimensionless i(x,Osb) function; AlGendy & Morsink et al 2014
jmom = imom*Osb*sqrt(c^2*R/G/M)
const wp = 2*jmom

# invariant quadrupole moment
qinv = quad + (4/3)*beta


println("beta: $beta | q: $quad | j: $jmom | qinv: $qinv")


#-------------------------------------------------- 

include("strig.jl") #spherical trigonometric functions and ellipsoid shape
include("bender.jl")
include("rtrace.jl")
include("beaming.jl")
include("radiation.jl")
include("plot2d.jl")
include("img.jl")
#include("cspot.jl")
