#Compute flux in polar coordinates

function radiation(Ir,
                   rad, chi,
                   phi, theta, cosa,
                   X, Xob, Osb, sini, earea)

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
    ####
    #vphi = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    #b = R*vphi/c
    #vw = Rgm*(1/enu)*sin(theta)*w #isoradial space vel
    #bp = R*vw/c
    #gamma = 1/sqrt(1 - b^2)
    #cosi = sqrt(1-sini^2)
    #sina = sqrt(1-cosa^2)
    #cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
    #cosz = -sina*sini*sin(phi)/sqrt(1-cospsi^2)
    #delta = (1/gamma)/(1 - b*cosz)
    #EEd = delta*enu*(1 + cosz*bp)
    #else
    
    vz = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    bz = R*vz/c

    gamma = 1/sqrt(1-bz^2)
    delta = 1/gamma/(1 + Lz*(2pi*fs)/(G*M/c^2))
    EEd = delta*enu
    #end
    
    dS = (Rgm)^2*sin(theta)*sqrt(1 + fa^2)
    cosap = cosa * delta
    dOmega = dS*cosap
    
    dF = (EEd^3)*Ir(cosap) * earea * delta
    #dF = (EEd^3)*dOmega*Ir(cosap)
    #dF = (EEd^3)*dOmega*delta
    
    return dF, EEd
    #return gamma, gamma2
    #return dS, dOmega
end

Areas = zeros(Nrad, Nchi)
Flux = zeros(Nrad, Nchi)
Reds = zeros(Nrad, Nchi)

#pts = [(0.0, 0.0),
#       (0.0, 0.0),
#       (0.0, 0.0),
#       (0.0, 0.0)]
pts_phis = zeros(4)
pts_thetas = zeros(4)
       
tic()
println("Computing radiation from the star")
for i = 2:Nchi-1
    chi = chi_grid[i]
    
    if mod(i,10) == 0; println("chi: $(round(chi/2pi,2))"); end
    
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
        if hits[j+1, i] < 1
            jni = [(j-1, i-1),
                   (j-1, i+1),
                   (j  , i)]
            iarea = abs(rad_grid[j] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
        else
            jni = [(j-1, i-1),
                   (j-1, i+1),
                   (j+1, i+1),
                   (j+1, i-1)]
            iarea = abs(rad_grid[j+1] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
        end
        
        #build corners
        q = 0
        for (j1, i1) in jni
            q += 1
            #pts[q] = (rad[j1,i1], chi_grid[j1,i1])
            pts_phis[q] = Phis[j1, i1]
            pts_thetas[q] = Thetas[j1, i1]
        end

        earea = area_sphere_lambert(phi, theta, pts_phis[1:q], pts_thetas[1:q], Rq, ecc)
        Areas[j,i] = earea/iarea
        
        dF, dE = radiation(Ir,
                           rad, chi,
                           phi, theta, cosa,
                           X, Xob, Osb, sini, earea/iarea)

        Flux[j,i] = dF
        Reds[j,i] = dE
                
    end#j over rad
end#i over chi

#Fill the middle point for r = 0 by taking the mean of surrounding values
meanflux = 0.0
meanreds = 0.0
meanarea = 0.0
meanN = 0
for i = 1:Nchi
    if Flux[2,i] >0.0
        meanflux += Flux[2,i]
        meanreds += Reds[2,i]
        meanarea += Areas[2,i]
        meanN += 1
    end
end

Flux[1,:] = meanflux/meanN
Reds[1,:] = meanreds/meanN
Areas[1,:]= meanarea/meanN

toc()

area_interp    = interpolate((rad_grid, chi_grid), Areas, method)
flux_interp    = interpolate((rad_grid, chi_grid), Flux, method)
reds_interp    = interpolate((rad_grid, chi_grid), Reds, method)
