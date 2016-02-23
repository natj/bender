using Winston
using Colors
using Interpolations

include("plot2d.jl")

######################
# Physical constants
const G = 6.67408e-8
const c = 2.99792458e10
const Msun = 1.98892e33
const km = 1.0e5
const ergkev = 6.2415e8 # erg/keV 
const cm_parsec = 3.2404e-23 #1cm/10kpc

#initial parameters in physical units
incl = deg2rad(60.0)
M    = 1.6Msun
R    = 15.0km
fs   = 600
#Dist = 1.0*cm_parsec


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
#const beta = 0.4454*Osb^2*X #Matter quadrupole moment; Morsink 2014
#const quad = -0.11*(Osb/X)^2 #Energy quadrupole moment; Morsink 2014
I = sqrt(X)*(1.136 - 2.53*X + 5.6*X^2) #Moment of inertia; Morsink 2014
#const wp = (2*I*(2pi*fs)/X^2) * ((0.728194*(M/Msun)^2)/(G*M/c^2))

const beta = 0.0
const quad = 0.0
const wp = 0.0
println("beta=$beta q=$quad wp=$wp")

#
wf = (2*I*2pi*fs/X^2)*X^3*(1-3*X)
println("w=$(wf/2pi)")


#Load spherical trigonometric functions and ellipsoid shape
include("strig.jl")



####################################
####################################
#Moments


function ptim(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)

    return exp(-2*x^3*nu2)*(2+x)^2*(-1 + a*x^3*(-1 + 3*x)*wp*sini)/(-2 + x)^2
end

#Radial moment p^r
function prad(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    rturn = false
    sq = (-16*(a^2 + b^2)*exp(4*x^3*nu2)*(-2 + x)^4*x^2)/(Rg^2*(2 + x)^4*(4 + (-1 + 4*B2)*x^2)^2) + (1 + a*x^3*(1 - x*3)*wp*sini^2)
    if sq < 0.0
        sq = -sq
        rturn = true
    end
    return -4e^(-zeta2*x^2)*sqrt(sq)/(x^2 - 4), rturn
end

#Theta moment p^θ
function pthe(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    tturn = false
    sq = a^2 + b^2 - a^2*(sini*csc(theta))^2
    if sq < 0.0
        sq = -sq
        tturn = true
    end
    
    return 16*exp(x^2*(2x*nu2-zeta2))*(x-2)*x^2*sqrt(sq)/(Rg^2*(2+x)^3*(4+(4*B2-1)*x^2)), tturn
end

#Phi moment p^ϕ
function pphi(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    if -1.0e-5 < theta < 1.0e-5
        return 0.0
    end

    return (16*x^2*((a*exp(4*x^3*nu2)*(-2 + x)^4*csc(theta)^2*sini)/(Rg^2*(4 + (-1 + 4*B2)*x^2)^2) - (x*(2 + x)^4*(1 - x*3)*wp*(1 + a*x^3*(1 - x*3)*wp*sini))/16.))/(exp(2*x^3*nu2)*(-4 + x^2)^2)
end

#Runge-Kutta step
function rk_step(rri, yni,
                 x, y, sini, wp, Rg)
    
    nu2   = beta/3.0 - quad*0.5*(3*cos(yni)^2-1)
    B2    = beta
    zeta2 = beta*(3*0.5*(3*cos(yni)^2-1)/4-1/3)

    pa        = ptim(x, y, sini, rri, nu2, B2, zeta2, wp, yni, Rg)
    pr, rturn = prad(x, y, sini, rri, nu2, B2, zeta2, wp, yni, Rg)
    pt, tturn = pthe(x, y, sini, rri, nu2, B2, zeta2, wp, yni, Rg)
    pp        = pphi(x, y, sini, rri, nu2, B2, zeta2, wp, yni, Rg)
    dr        = -Rg/rri^2 #dx to dr

    dt = dr*pa/pr
    dy = dr*pt/pr
    dz = dr*pp/pr
    
    return dt, dy, dz, tturn, rturn
end

#polar coordinate wrapper
function bender3p(rad, chi, sini,
                  X, Osb,
                  beta, quad, wp, Rg)

    #x = rad*sin(chi)
    #y = rad*cos(chi)
    #println("x=$x y=$y")
    
    bender3(rad*sin(chi), rad*cos(chi), sini,
            X, Osb,
            beta, quad, wp, Rg)
end

######################
function bender3(x, y, sini,
                 X, Osb,
                 beta, quad, wp, Rg)
    
    if y >= 0.0
        psign = -1.0
    else
        psign = 1.0
    end
    rsign = 1.0

    ######################
    #leapfrog tracing
    const h = 0.002
    const tol = 5.0e-7 #relative tol
    
    hit = true
    rr = 0.0

    #first step
    rr += h/256
    
    tm1 = 0.0
    ym1 = incl
    zm1 = 0.0

    #series expanded initial values for moments assuming O(rr)^1
    tn = tm1 + rr*(x^2 + y^2)/2/Rg #referenced from x=y=0 photon
    yn = ym1 - rr*y*psign
    zn = zm1 + rr*(x/sini/Rg)

    #initial values for integration
    hi = h
    level = 2.0^7

    #Photon trace arrays
    tns = Float64[tm1,tn]
    yns = Float64[ym1,yn]
    zns = Float64[zm1,zn]
    rns = Float64[0.0,rr]
    ers = Float64[0.0,0.0]
    lvs = Float64[1.0,1.0]

    Xob = 100.0
    maxr = rr
    
    oneturn = true
    while rr <= Xob
        err = 1.0
        
        yni = yn
        psigni = psign
        rsigni = rsign

        
        
        tp1 = 0.0
        yp1 = 0.0
        zp1 = 0.0

        ###XXX: debug variables
        #k1y = 0.0
        #k2y = 0.0
        #k1z = 0.0
        #k2z = 0.0
        #k1t = 0.0
        #k2t = 0.0
                    
        while err > tol && level <= 128.0
            level *= 2.0
            #level = max(level*2.0, min(128.0, 1.1*err/tol))
                        
            #reset current state
            yn = yni
            psign = psigni
            rsign = rsigni
            hi = h / level
            
            #take a step
            k1t, k1y, k1z, tturn1, rturn1 = rk_step(rr, yn, x, y, sini, wp, Rg)
            k2t, k2y, k2z, tturn2, rturn2 = rk_step(rr+rsign*hi, yn+ psign*k1y*hi, x, y, sini, wp, Rg)
            
            #check if our photon turns around
            #if (tturn1 || tturn2)
            #    psign *= -1.0
            #end
            #if (rturn1 || rturn2)
            #    rsign *= -1.0
            #end
            if tturn1
                psign *= -1.0
            end
            if rturn1
                rsign *= -1.0
            end

            
            #if rturn1 && rturn2 #|| rturn3 || rturn4 || rturn5 || rturn6
            #    hit = false
            #    #return rns, yns, zns, ers, lvs
            #    return tp1, yp1, zp1, Xob, hit, 0.0
            #end

            #rk21 adaptive
            tp1_o = tn + hi*k1t
            yp1_o = yn + hi*k1y*psign
            zp1_o = zn + hi*k1z

            tp1 = tn + hi*(0.5*k1t + 0.5*k2t)
            yp1 = yn + hi*(0.5*k1y + 0.5*k2y)*psign
            zp1 = zn + hi*(0.5*k1z + 0.5*k2z)

            #Error
            errt = abs(tp1 - tp1_o)
            erry = abs(yp1 - yp1_o)
            errz = abs(zp1 - zp1_o)

            err = max(abs(erry/yp1), abs(errz/zp1), 10*abs(errt/tp1)) #rel err
            #err = max(abs(erry), abs(errz)) #abs err
        end

        rr = rr + rsign*hi
        
        push!(lvs, log2(level))
        #push!(lvs, level*psign)
        push!(ers, err)

        level = max(0.5, level/4.0)
        
        
        #theta & theta
        tn = tp1
        yn = yp1
        zn = zp1

        #store photon trace
        push!(tns, tn)
        push!(yns, yn)
        push!(rns, rr)
        push!(zns, zn)
        #push!(zns, k1z)
        

        nu2   = beta/3.0 - quad*0.5*(3*cos(yn)^2-1)
        B2    = beta
        enu = (1-rr/2)/(1+rr/2)*exp(nu2*rr^3)
        B = (1-rr/2)*(1+rr/2) + B2*rr^2

        #Radius (isoradial)
        Rgm, dtR = Rgmf(yn, X, Osb)
        Xobi = X/Rgm
        Xob = Xobi*B/enu #isotropic x to be referenced with rr

        #Radius (isotropic)
        #Rgm, dtR = Rgmf(yn, X, Osb)
        #Xobi = X/Rgm
        #Xob = Xobi #*enu/B #isotropic x; assuming spherical star, i.e. no conversion

        
        #Keep track of photon U-turns
        if rr > maxr
            maxr = rr
        end

        #Terminate photon if we miss by 5%
        if rr < maxr
            if rr < Xob/0.95
                hit = false
                break
            end
        end

        #Break down if we are close to the star;
        #helps with surface detection and with 1/r^3 functions
        if rr > Xob*0.95
            level = 128.0
        end
    end

    time = tn
    theta = yn
    phi = mod2pi(pi-zn)-pi

    #Emission regions (in comparison to Cadeau 2007
    #return time, theta, phi, Xob, hit, 0.0
    
    #psignpm = psign*sign(cos(theta))
    #time = 0
    #if rsign < 0 && psignpm < 0
    #    time = 1
    #elseif rsign > 0 && psignpm < 0
    #    time = 2
    #elseif rsign > 0 && psignpm > 0
    #    time = 3
    #else
    #elseif rsign > 0 && psignpm > 0
    #    time = 4
    #end

    #return time, psign, rsign, Xob, hit, 0.0
    
    if !hit
        return time, theta, phi, Xob, hit, 0.0
    end

    
    #Emission angle alpha
    ##########
    nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
    B2    = beta
    zeta2 = beta*(3*0.5*(3*cos(theta)^2-1)/4-1/3)
    Rgm, dR = Rgmf(theta, X, Osb)
    
    enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
    B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
    ezeta = (1-Xob/2)*(1+Xob/2) + zeta2*Xob^2
    C = (x^2 + y^2)
    Lz = x*sini
    w = wp*Xob^3*(1-3*Xob)
    
    fa = (B/enu/ezeta)*dR/Rgm
    cosg = 1/sqrt(1 + fa^2)
    sing = fa*cosg

    #cosg = 0.5
    #sing = 0.5
    
    #Cadeau 2007 hit test
    #if x > 0
    #    pm = -1
    #else
    #    pm = 1
    #end
    #cad = (1/(B*enu^2)*sin(theta))/Xob/sini/(1 + pm*w/(B*enu^2)*sin(theta)/Xob)
    #cad = (1/(B*enu^2)*sin(theta))/Xob/sini
    #if abs(x) > cad
    #    hit = false
    #    return time, theta, phi, Xob, hit, 2.0
    #end

    
    ####
    ep = 1.0
    ds = 1.0
    mult = ds*(1 - (B*w*sin(theta)*Rg/Xob/enu^2)^2)/ep

    t1 = (ep + Lz*w)
    
    sq2 =  (ep + Lz*w)^2 - C*(enu^4)*(Xob^2)/(B^2)/(Rg^2)
    if sq2 < 0; sq2 = 0.0; end
    t2 = cosg*sqrt(sq2)
    
    sq3 = C - Lz^2*csc(theta)^2
    if sq3 < 0; sq3 = 0.0; end
    t3 = psign*sing*enu^2/B*(Xob/Rg)*sqrt(sq3)
    
    cosa = 1 + mult*(-t1 + t2 + t3)
    
    if cosa < 0
        hit = false
    end
    cosa = clamp(cosa, 0.0, 1.0)
    
    #return rns, yns, zns, ers, lvs, hit
    return time, phi, theta, Xob, hit, cosa
end
