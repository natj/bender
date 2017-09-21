

####################################
####################################
#Moments

#time moment p^t
function ptim(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)

    # println("enu0:",(1-x/2)/(1+x/2))
    # println("enu2:",exp(nu2*x^3))
    # println("enu:",((1-x/2)/(1+x/2)*exp(nu2*x^3))^2)
    # println("w:",x^3*(1 - 3*x)*wp)
    # println(" enu2:",exp(-2*x^3*nu2) * (2 + x)^2 / (-2 + x)^2)
    # println(" w   :",(-1 + a*x^3*(-1 + 3*x)*wp*sini))


    return exp(-2*x^3*nu2) * (2 + x)^2 * (-1 + a*x^3*(-1 + 3*x)*wp*sini)/(-2 + x)^2
end

#Radial moment p^r
function prad(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    rturn = false
    sq = (-16*(a^2 + b^2)*exp(4*x^3*nu2)*(-2 + x)^4*x^2)/(Rg^2*(2 + x)^4*(4 + (-1 + 4*B2)*x^2)^2) + (-1 + a*x^3*(-1 + x*3)*wp*sini)^2
    #sq1 = (-16*(a^2 + b^2)*exp(4*x^3*nu2)*(-2 + x)^4*x^2)/(Rg^2*(2 + x)^4*(4 + (-1 + 4*B2)*x^2)^2)
    #sq2 = (1 - a*x^3*(-1 + x*3)*wp*sini)^2
    #println("rad sq ",sq)
    #println("rad sq1 ",sq1)
    #println("rad sq2 ",sq2)
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
    #println("the sq",sq)
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
function rk_step(rrs, yis,
                 x, y, sini, wp, Rg)
    
    nu2   = beta/3.0 - quad*0.5*(3*cos(yis)^2-1)
    B2    = beta
    zeta2 = beta*(3*0.5*(3*cos(yis)^2-1)/4-1/3)

    pa        = ptim(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
    pr, rturn = prad(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
    pt, tturn = pthe(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
    pp        = pphi(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
    dr        = -Rg/rrs^2 #dx to dr

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

    #println(x," ",y)
    
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

    #
    tni = 0.0
    yni = 0.0
    zni = 0.0
    rri = 0.0
    xoi = 0.0
    maxri = 0.0
    rsigni = rsign
    psigni = psign

    counter = 0
    
    Xob = 100.0
    maxr = rr
    
    #println()
    #println("tn ", tn)
    #println("yn ", yn)
    #println("zn ", zn)
    #println()



    oneturn = true
    #while rr <= Xob
    while true
        #println(rr)
        
        err = 1.0

        tni = tn
        yni = yn
        zni = zn
        rri = rr
        xoi = Xob
        maxri = maxr
        
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

        level_break = true
        #while err > tol #&& level <= 256.0
        while err > tol && level_break
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

            #println("tp1 ",tp1)
            #println("yp1 ",yp1)
            #println("zp1 ",zp1)

            if level > 512.0
                level_break = false
            end
        end

        rr += rsign*hi
        
        #push!(lvs, log2(level)*rsign)
        #push!(lvs, level*psign)
        #push!(ers, err)

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
        ##push!(zns, k1z)

        nu2   = beta/3.0 - quad*0.5*(3*cos(yn)^2-1)
        B2    = beta
        #nu2   = 0.0
        #B2    = 0.0
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

        #println("$tn, $yn, $zn $err $Rgm $Xob")
        
        #Keep track of photon U-turns
        if rr > maxr
            maxr = rr
        end

        #Terminate photon if we miss by 5%
        if rr < maxr
            if rr < Xob/0.95
                #println("terminate")
                hit = false
                break
            end
        end

        #Slow down if we are close to the star;
        #helps with surface detection and with 1/r^3 functions
        #if rr > Xob*0.95
        #    level = 128.0
        #end

        #serr = (Xob - rr)/Xob
        #if rr > 0.18
        #    println("Xob: ",Xob, " rr:", rr, " serr: ",serr, " lvl:", level, " hit: ",hit, " yni:", yn)
        #end

        #if false
        if rr >= Xob
            #println("    penetration")
            
            serr = (Xob - rr)/Xob
            counter += 1

            if counter > 50
                break
            end
            
            if abs(serr) > 1.0e-6
                #println("set previous step")
                
                #previous (non-penetrated) step
                tn = tni
                yn = yni
                zn = zni
                rr = rri
                maxr = rri
                #maxr = maxri
                #Xob = xoi
                rsign = rsigni
                psign = psigni
                level *= 4.0

                #level = min(level, 128.0)
            else
                break
            end
        end
    end

    #Henon's trick
    
    #println("Xob: ",Xob)
    #println("rr: ",rr)
    
    err = (Xob - rr)/Xob
    stol = 1.0e-8 #surface detection tolerance
    #if !hit
    #    println("aaaaaaaaaaaaa")
    #end
    #if hit
    if false

        #set boundaries
        rmin = rri
        rdmin = xoi - rmin
        
        rmax = rr
        rdmax = Xob - rmax

        println("l: ",length(yns))
        println("rmin: ",rmin)
        println("rdmin: ",rdmin)
        println("rmax: ",rmax)
        println("rdmax: ",rdmax)
        
        #previous (non-penetrated) step
        #tn = tni
        #yn = yni
        #zn = zni
        #rr = rri
        #Xob = xoi
        #rsign = rsigni
        #psign = psigni
        
        tp1 = 0.0
        yp1 = 0.0
        zp1 = 0.0
        rro = rr
        
        while abs(err) > stol
        #for j = 0:20
            
            hi = Xob - rr
            #hi = clamp(hi, -h/64, h/64)
            
            #take a step
            k1t, k1y, k1z, tturn1, rturn1 = rk_step(rr, yn, x, y, sini, wp, Rg)
            k2t, k2y, k2z, tturn2, rturn2 = rk_step(rr+rsign*hi, yn+ psign*k1y*hi, x, y, sini, wp, Rg)

            if tturn1; psign *= -1.0; end
            if rturn1; rsign *= -1.0; end

            #rk1 (euler)
            #tp1 = tn + hi*k1t
            #yp1 = yn + hi*k1y*psign
            #zp1 = zn + hi*k1z
            
            #rk21 adaptive
            tp1 = tn + hi*(0.5*k1t + 0.5*k2t)
            yp1 = yn + hi*(0.5*k1y + 0.5*k2y)*psign
            zp1 = zn + hi*(0.5*k1z + 0.5*k2z)

            rro = rr + rsign*hi

            #compute new Xob
            nu2   = beta/3.0 - quad*0.5*(3*cos(yp1)^2-1)
            B2    = beta
            enu = (1-rro/2)/(1+rro/2)*exp(nu2*rro^3)
            B = (1-rro/2)*(1+rro/2) + B2*rro^2

            Rgm, dtR = Rgmf(yp1, X, Osb) #radius (isoradial)
            Xobi = X/Rgm
            Xob = Xobi*B/enu #isotropic x to be referenced with rro

            err = (Xob - rro)/Xob
            #println("err: ", err, " hi: ",hi)

            
            #update step
            #rr = rro
            #tn = tp1
            #yn = yp1
            #zn = zp1

           
            #store photon trace
            #push!(tns, tp1)
            #push!(yns, yp1)
            #push!(rns, rro)
            #push!(zns, zp1)
            #push!(ers, err)
        end

        #push!(tns, tp1)
        #push!(yns, yp1)
        #push!(rns, rro)
        #push!(zns, zp1)
        #push!(ers, err)

        tn = tp1
        yn = yp1
        zn = zp1
        #rr = rr
    end

    #println()
    #println("l2: ",length(yns))
    #println("Xob2: ",Xob)
    #println("rr2: ",rr)
    #println()
    
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
    
    #if !hit
    #    return time, theta, phi, Xob, hit, 0.0
    #end

    
    #Emission angle alpha
    ##########
    nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
    B2    = beta
    zeta2 = beta*((4/3)*0.5*(3*cos(theta)^2 - 1) - 1/3)
    Rgm, dR = Rgmf(theta, X, Osb)
    
    enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
    B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
    ezeta = (1-Xob/2)*(1+Xob/2)*exp(zeta2*Xob^2)
    C = (x^2 + y^2)
    Lz = x*sini
    w = wp*Xob^3*(1-3*Xob)
    
    fa = (B/enu/ezeta)*dR/Rgm
    cosg = 1/sqrt(1 + fa^2)
    sing = fa*cosg

    vz = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    bz = R*vz/c
    gamma = 1/sqrt(1 - bz^2)
    eta =  1/(1 + Lz*(2pi*fs)*(G*M/c^3))
    delta = (eta/gamma)
    EEd = delta*enu

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
    cosao = 1 + mult*(-t1 + t2 + t3)
    

    #new emission angle using projection operator definition
    # pr 
    sq2 =  (ep + Lz*w)^2 - C*(enu^4)*(Xob^2)/(B^2)/(Rg^2)
    if sq2 < 0; sq2 = 0.0; end
    pr = sqrt(sq2)*ezeta/enu^2


    sq3 = C - Lz^2*csc(theta)^2
    if sq3 < 0; sq3 = 0.0; end
    pt = psign*sqrt(sq3)*ezeta/B

    dotpr = (cosg*pr + sing*pt*Xob)
    cosa = enu^2/ezeta*delta*dotpr

    #println("cosao: $cosao | cosa: $cosa | ratio $(cosao/cosa)")

    if cosa < 0
        hit = false
    end
    cosa = clamp(cosa, 0.0, 1.0)
    
    #return rns, yns, zns, tns, ers, lvs, hit
    return time, phi, theta, Xob, hit, cosa
end
