#Compute flux in polar coordinates
include("bender.jl")


#Planck function B_E(T) = erg cm^-2 s^-1 str^-1 keV^-1
#constbb = hzkeV^4 * 2*h/c^2
#hzkeV = 2.41789894e17
BE(T, Ener) = 5.039617e22 * (Ener .^3.0)./expm1(Ener ./ T)

#Planck function photon flux N_E(T) = ph cm^-2 s^-1 str^-1 keV^-1
#constbb = hzkeV^4 * 2*h/c^2
#hzkeV = 2.41789894e17
NBE(T, Ener) = 5.039617e22 * (Ener .^3.0)./expm1(Ener ./ T) ./ Ener * ergkev

#Number of black body photons N(T) = \int dE B_E(T)/E(keV) * ergkev
#\int dE E^2 /(exp(E/T)-1) = Gamma(3)Zeta(3) T^3 
NB(T) = 5.039617e22*6.2415e8*2.404*T.^3.0

#Energy of black body photons E(T) = \int dE B_E(T)
#\int dE E^3 /(exp(E/T)-1) = Gamma(4)Zeta(4) T^4 = 6 * pi^4/90  * T^4
EB(T) = 5.039617e22*(pi^4/15.0)*T.^4 * ergkev

#Beaming function
Beam(mu) = 1.0 #Lambertian
#Beam(mu) = 0.42822 + 0.92236*mu - 0.085751*mu^2 #approx Hopf


#Compute blackbody flux elements
function bbfluxes(EEd, delta, cosa)

    dist = 1.0 / cm_parsec #10 kpc
    d2 = dist^2
    
    const Energs = [2.0, 6.0, 12.0]
    const Teff = 2.0
    
    #Collect flux for different energies
    fluxE = (EEd)^3 .* BE(Teff, Energs ./ EEd) .* Beam(cosa*delta) *delta #energy flux
    fluxNE = (EEd)^2 .* NBE(Teff, Energs ./ EEd) .* Beam(cosa*delta) *delta #photon flux
    
    #Bolometric fluxes
    fluxNB = (EEd)^3 * NB(Teff) * Beam(cosa*delta) *delta #photon bol flux
    fluxB = (EEd)^4 * EB(Teff) * Beam(cosa*delta) *delta #energy bol flux

    return fluxE ./ d2, fluxNE ./ d2, fluxNB/d2, fluxB/d2
end



function radiation(rad, chi,
                   phi, theta, cosa,
                   X, Xob, Osb, sini)

    nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
    B2    = beta
    zeta2 = beta*(3*0.5*(3*cos(theta)^2-1)/4-1/3)
    Rgm, dR = Rgmf(theta, X, Osb)

    enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
    B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
    ezeta = (1-Xob/2)*(1+Xob/2) + zeta2*Xob^2

    C = rad^2
    Lz = sini*sqrt(C)*sin(chi)
    w = wp*Xob^3*(1-3*Xob)

    fa = (B/enu/ezeta)*dR/Rgm
    
    cosg = 1/sqrt(1 + fa^2)
    sing = fa*cosg


    #if false
    #########################
    vphi = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    b = R*vphi/c
    vw = Rgm*(1/enu)*sin(theta)*w #isoradial space vel
    bp = R*vw/c
    gamma = 1/sqrt(1 - b^2)
    cosi = sqrt(1-sini^2)
    sina = sqrt(1-cosa^2)
    cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
    cosz = -sina*sini*sin(phi)/sqrt(1-cospsi^2)
    delta = (1/gamma)/(1 - b*cosz)
    EEd = delta*enu #*(1 + cosz*bp)
    delta2 = delta
    #else


    #########################
    ##vz = Rgm*(1/enu)*sin(theta)*(2pi*fs) #isoradial velo (wrong)
    vz = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo 
    bz = R*vz/c

    gamma = 1/sqrt(1 - bz^2)
    dtaudt = (enu^2)/gamma
    delta = (1/gamma)/(1 + Lz*(2pi*fs)/(G*M/c^2)/dtaudt)
    EEd = delta*enu


    #end
    
    #dS = (Rgm)^2*sin(theta)*sqrt(1 + fa^2)
    #cosap = cosa * delta #abberration
    #dOmega = dS*cosap

    #dF = (EEd^4)*Ir(cosap)*earea #bolometric
    #dF = (EEd^3)*Ir(cosap)*earea #monochromatic
    #dF = (EEd^3)*dOmega*Ir(cosap)
    #dF = (EEd^3)*dOmega*delta

    #return delta2, delta
    return EEd, delta

    #return dF, EEd
    #return gamma, gamma2
    #return dS, dOmega
end

#Areas = zeros(Nrad, Nchi)
#dFlux = zeros(Nrad, Nchi)
Reds = zeros(Nrad, Nchi)
Deltas = zeros(Nrad, Nchi)

#Additional F_E and F
#NE = 3
#Energs = [2.0, 6.0, 12.0]

#FluxE = zeros(Nrad, Nchi, NE)
#FluxB = zeros(Nrad, Nchi)


#
#pts = [(0.0, 0.0),
#       (0.0, 0.0),
#       (0.0, 0.0),
#       (0.0, 0.0)]
pts_phis = zeros(4)
pts_thetas = zeros(4)
       
tic()
print("Computing radiation from the star...")
for i = 2:Nchi-1
    chi = chi_grid[i]
    
    #if mod(i,10) == 0; println("chi: $(round(chi/2pi,2))"); end
    
    for j = 2:Nrad-1
        rad = rad_grid[j]

        if hits[j, i] < 1; break; end

    
        #Ray traced photons
        ####
        phi = Phis[j, i]
        theta = Thetas[j, i]
        time = Times[j, i]
        Xob = Xs[j, i]
        cosa = cosas[j, i]

        #emitting area size
        ###

        #define corners depending on if we are on the edge
        #if hits[j+1, i] < 1
        #    jni = [(j-1, i-1),
        #           (j-1, i+1),
        #           (j  , i)]
        #    iarea = abs(rad_grid[j] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
        #else
        #    jni = [(j-1, i-1),
        #           (j-1, i+1),
        #           (j+1, i+1),
        #           (j+1, i-1)]
        #    iarea = abs(rad_grid[j+1] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
        #end
        
        #build corners
        #q = 0
        #for (j1, i1) in jni
        #    q += 1
        #    #pts[q] = (rad[j1,i1], chi_grid[j1,i1])
        #    pts_phis[q] = Phis[j1, i1]
        #    pts_thetas[q] = Thetas[j1, i1]
        #end
        #sarea = area_sphere_lambert(phi, theta, pts_phis[1:q], pts_thetas[1:q], Rq, ecc)
        #Areas[j,i] = sarea/iarea
        #Areas[j,i] = iarea
        #Areas[j, i] = 1.0
        
        #Ir(cosa) = 1.0 #lambertian intensity
        
        EEd, delta = radiation(rad, chi,
                               phi, theta, cosa,
                               X, Xob, Osb, sini)

        
        #dFlux[j, i] = dF
        Reds[j, i] = EEd
        Deltas[j, i] = delta
                
    end#j over rad
end#i over chi

#Fill the middle point for r = 0 by taking the mean of surrounding values
#meanflux = 0.0
meanreds = 0.0
meandelta = 0.0
#meanarea = 0.0
meanN = 0
for i = 1:Nchi
    if Reds[2,i] >0.0
        #meanflux += dFlux[2,i]
        meanreds += Reds[2,i]
        meandelta += Deltas[2,i]
        #meanarea += Areas[2,i]
        meanN += 1
    end
end

#dFlux[1,:] = meanflux/meanN
Reds[1,:] = meanreds/meanN
Deltas[1,:] = meandelta/meanN
#Areas[1,:]= meanarea/meanN

toc()

#area_interp    = interpolate((rad_grid, chi_grid), Areas, method)
#flux_interp    = interpolate((rad_grid, chi_grid), dFlux, method)
reds_interp    = interpolate((rad_grid, chi_grid), Reds, method)
delta_interp    = interpolate((rad_grid, chi_grid), Deltas, method)
