using Winston

include("plot2d.jl")


const G = 6.67384e-8
const c = 2.99792458e10
const Msun = 1.9885469e33 #XXX
const km = 1.0e5
const ergkev = 6.24150934326e8 # erg/keV 
const cm_parsec = 3.24077929e-23 #1cm/10kpc # 3.08567758135
const constbb = 5.040366e22



fs   = 700

rgrid = linspace(8.0, 18, 100)
mgrid = linspace(1.0, 2.0, 100)


qgrid = zeros(length(mgrid), length(rgrid))

function quad_F(R,M,fs)
		X = G*M/(R*c^2)
                Osb = (2pi*fs)*sqrt(R^3/(G*M))
		quad = -0.11*(Osb/X)^2
	return quad
end

for xr = 1:length(rgrid)
        R = rgrid[xr]*km
	for yr = 1:length(mgrid) 
		M = mgrid[yr]*Msun

		Xob = G*M/(R*c^2)
                Osb = (2pi*fs)*sqrt(R^3/(G*M))
		#quad = -0.11*(Osb/X)^2
		X = G*M/(R*c^2)
		beta = 0.4454*5*Osb^2*X #Matter quadrupole moment; AlGendy & Morsink 2014

		quad = quad_F(R,M,fs)
                theta = 048
                nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
                enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)

		val = enu/((1-Xob/2)/(1+Xob/2))

                qgrid[yr, xr] = val
	end
end


p1 = plot2d(qgrid, collect(rgrid),
 		   collect(mgrid),
		   0,
                   0.0,
		   0.0,
                   "Blues",
                   xlabel="Radius", 
                   ylabel="Mass")
