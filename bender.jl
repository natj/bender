using Winston
using Colors
using Grid

include("plot2d.jl")

######################
# Physical constants
const G = 6.67259e-8
const c = 2.99792458e10
const Msun = 1.99e33
const km = 1.0e5

#initial parameters in physical units
incl = pi/4
#M    = 1.65Msun #5.0
M    = 1.4Msun
R    = 10km
fs   = 50

#initial step with flat space
#increase level around ptheta zero

#derived dimensionless values
const sini = sin(incl)
#const Rg = G*M/c^2
const Rg = 1.0
#const Rg = 1.4774272698670399
#const Rg = M*(G/c^2)
const X = G*M/(R*c^2)
const Osb = (2pi*fs)*sqrt(R^3/(G*M))
const U = 2*G*M/(R*c^2)
println("x=$X ($U) and Osb=$Osb incl=$incl")

#Hartle-Thorne parameters
#const beta = 0.4454*Osb^2*X #Matter quadrupole moment; Morsink 2014
#const quad = -0.11*(Osb/X)^2 #Energy quadrupole moment; Morsink 2014
I = sqrt(X)*(1.136-2.53*X+5.6*X^2) #Moment of inertia; Morsink 2014
const wp = 2*I*(2pi*fs)/X^2 / (G*M/c^2)

const beta = 0.0
const quad = 0.0
#const wp = 0.0
println("beta=$beta q=$quad wp=$wp")


#
wf = (2*I*2pi*fs/X^2)*X^3*(1-3*X)
println("w=$(wf/2pi)")

######################
######################
#Radial function by AlGendy & Morsink 2014
#Circumferential radius
function Rgmf(theta, X, Osb)
    const o20 = -0.788
    const o21 = 1.030

    #Radial function
    req=1.0
    o2=(o20+o21*X)*Osb^2
    Rgm = req*(1+o2*cos(theta)^2)

    #derivative dR/dtheta
    dtR = -2*Osb^2*(o20+o21*X)*cos(theta)*sin(theta) #from AlGendy & Morsink 2014

    #return 1.0, 0.0
    return Rgm, dtR
end

###################
# Trigonometric functions

eRe, dR = Rgmf(pi/2, X, Osb) #equatorial radius
eRp, dR = Rgmf(0, X, Osb) #pole radius
const ecc = sqrt(1 - (eRp/eRe)^2) #eccentricity
if ecc != 0.0
    const qp = 1 + (1-ecc^2)/(2*ecc)*log((1+ecc)/(1-ecc)) #authalic pole latitude
    const Rq = eRe*sqrt(qp/2) #authalic radius
else
    const qp = 0.0
    const Rq = 1.0
end
println("ecc = $ecc Rq = $Rq")

#authalic colatitude 
function authalic_lat(colat, ecc, qp)

    lat = pi/2 - colat
    sinl = sin(lat)
    q = (1-ecc^2)*(sinl/(1 - (ecc*sinl)^2) - log((1 - ecc*sinl)/(1 + ecc*sinl))/2/ecc)
    alat = asin(q/qp)
    acolat = pi/2 - alat
    
    return acolat
end

#great circle distance (on a sphere)
function great_circle_dist(lon1, lon2, col1, col2)

    lat1 = pi/2 - col1
    lat2 = pi/2 - col2
    dlon = lon2 - lon1
    
    xx = (cos(lat1)*sin(dlon))^2 + (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))^2
    yy = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)

    return atan2(sqrt(xx), yy)
end

#Area of a polygon on a sphere/spheroid
function area_sphere(phi0,the0,lons, cols, Rarea=1, ecc=0.0)

    N = length(cols)
    
    #transform to authalic sphere if we have spheroid
    if ecc != 0.0
        for i = 1:N
            cols[i] = authalic_lat(cols[i], ecc, qp)
        end
        the0 = authalic_lat(the0, ecc, qp)
    end

    #sigs = zeros(N)
    #for i = 1:N-1
    #    sigs[i] = great_circle_dist(lons[i], lons[i+1], cols[i], cols[i+1]) 
    #end
    #sigs[N] = great_circle_dist(lons[N], lons[1], cols[N], cols[1]) 
    #area = Rarea^2 * (sum(abs(sigs)) - (N-2)*pi)

    #area = (lons[2] - lons[N])*cos(cols[1])
    #for i = 2:N-1
    #    area += (lons[i+1] - lons[i-1])*cos(cols[i])
    #end
    #area += (lons[1] - lons[N-1])*cos(cols[N])
    #area = (area*R^2)/2

    
    function azimuth(lon1, lon2, col1, col2)
        lat1 = col1 - pi/2
        lat2 = col2 - pi/2
        
        f1 = cos(lat2) * sin(lon2-lon1)
        f2 = cos(lat1) * sin(lat2)
        f3 = sin(lat1) * cos(lat2) * cos(lon2-lon1)
        az = atan2(f1, f2-f3)
        return mod2pi(az)
    end

    #cols = pi/2 - cols
    #lats = cols - pi/2
    #lat0 = 0

    col0 = pi/2 - the0 
    lon0 = phi0
    #col0 = pi/2
    #lon0 = 0

    colat = zeros(N+1)
    azimu = zeros(N+1)
    for i = 1:N
        colat[i] = great_circle_dist(lon0, lons[i], col0, cols[i])
        azimu[i] = azimuth(lon0, lons[i], col0, cols[i])
    end
    colat[N+1] = great_circle_dist(lon0, lons[1], col0, cols[1])
    azimu[N+1] = azimuth(lon0, lons[1], col0, cols[1])
    
    #step size
    daz = diff(azimu)
    daz = pi .* ((abs(daz) ./ pi) - 2.* ceil(((abs(daz) ./ pi) .- 1) ./ 2 )) .* sign(daz) #wrap to [-pi,pi]

    #average surface distance
    deltas = diff(colat) ./ 2
    colat = colat[1:N] .+ deltas

    integ = (1 .- cos(colat)) .* daz
    area = abs(sum(integ))/4/pi
    area = min(area, 1-area)

    return area*4*pi*Rarea
end

function area_sphere_tri(lons, cols, Rarea=1, ecc=0.0)
    N = length(cols)
    
    #transform to authalic sphere if we have spheroid
    if ecc != 0.0
        for i = 1:N
            cols[i] = authalic_lat(cols[i], ecc, qp)
        end
    end

    sigs = zeros(N)
    for i = 1:N-1
        sigs[i] = great_circle_dist(lons[i], lons[i+1], cols[i], cols[i+1])
    end
    sigs[N] = great_circle_dist(lons[N], lons[1], cols[N], cols[1]) 
    area = Rarea^2 * (sum(abs(sigs)) - (N-2)*pi)

    return area
end



######
#eqt = pi/2
#nu2   = beta/3.0 - quad*0.5*(3*cos(eqt)^2-1)
#B2    = beta
#zeta2 = beta*(3*0.5*(3*cos(eqt)^2-1)/4-1/3)

#Rgm, dR = Rgmf(eqt, X, Osb)
#Xob = X/Rgm

#enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
#B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
#isoR = (B/enu)*R
#const U = 2*G*M/(isoR*c^2)

#println("R=$(round(R/km,2))")
#println("isoR = $(round(isoR/km,2)) U=$U")

#Moments

function ptim(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)

    return exp(-2*x^3*nu2)*(2+x)^2*(-1 + a*x^3*(-1 + 3*x)*wp*sini)/(-2 + x)^2
end

#Radial moment pr
function prad(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    rturn = false
    sq = (-16*(a^2 + b^2)*exp(4*x^3*nu2)*(-2 + x)^4*x^2)/(Rg^2*(2 + x)^4*(4 + (-1 + 4*B2)*x^2)^2) + (1 + a*x^3*(1 - x*3)*wp*sini^2)
    if sq < 0.0
        sq = -sq
        #sq = 0.001
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
        #sq = 0.0
        tturn = true
    end

    return 16*exp(x^2*(2x*nu2-zeta2))*(x-2)*x^2*sqrt(sq)/(Rg^2*(2+x)^3*(4+(4*B2-1)*x^2)), tturn
end

#Phi moment p^ϕ
function pphi(a, b, sini,
              x, nu2, B2, zeta2, wp, theta, Rg)
    
    if -1e-5 < theta < 1e-5
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


######################
function bender3(x, y, sini,
                 X, Osb,
                 beta, quad, wp, Rg)
    
    if y >= 0.0
        psign = -1.0 # XXX: TODO better psign initial val
    else
        psign = 1.0
    end
    rsign = 1.0
    
    ######################
    #leapfrog tracing
    const h = 0.002
    #const tol = 0.0003 #absolute tol
    #const tol = 1.0e-7 #relative tol
    const tol = 1.0e-5 #relative tol
    
    hit = true
    rr = 0.0

    #first step
    rr += h/256
    
    tm1 = 0.0
    ym1 = incl
    zm1 = 0.0
    dt, dy, dz, tturn1, rturn1 = rk_step(rr, ym1, x, y, sini, wp, Rg)

    tn = tm1 + h*dt
    yn = ym1 + h*dy*psign
    zn = zm1 + h*dz

    #initial values for integration
    #rr -= h/256
    hi = h
    level = 0.5

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
        #rr += hi
        #rr = rr + rsign*hi
        
        err = 1.0
        
        yni = yn
        psigni = psign
        rsigni = rsign

        
        
        tp1 = 0.0
        yp1 = 0.0
        zp1 = 0.0

        ###XXX
        #k1y = 0.0
        #k2y = 0.0
        #k1z = 0.0
        #k2z = 0.0
        #k1t = 0.0
        #k2t = 0.0
            
        
        while err > tol && level <= 2^7
            level *= 2
            
            #reset current state
            yn = yni
            psign = psigni
            rsign = rsigni
            hi = h / level
            
            #take a step
            k1t, k1y, k1z, tturn1, rturn1 = rk_step(rr, yn, x, y, sini, wp, Rg)
            k2t, k2y, k2z, tturn2, rturn2 = rk_step(rr+rsign*hi, yn+ psign*k1y*hi, x, y, sini, wp, Rg)
            
            #check if our photon turns around
            if (tturn1 || tturn2)
                psign *= -1.0
            end
            if (rturn1 || rturn2)
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
            #errt = 0.0
            erry = abs(yp1 - yp1_o)
            errz = abs(zp1 - zp1_o)

            err = max(abs(erry/yp1), abs(errz/zp1), abs(errt/tp1)) #rel err
            #err = max(abs(erry), abs(errz)) #abs err
            #level *= 2
        end

        rr = rr + rsign*hi
        
        push!(lvs, log2(level))
        #push!(lvs, level*psign)
        push!(ers, err)

        
        #if level >= 5
        #    level = level/2.0
        #else
        #    level = 1.0
        #end
        #level = level/2.0
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
        #push!(zns, k2y)
        

        nu2   = beta/3.0 - quad*0.5*(3*cos(yn)^2-1)
        B2    = beta
        enu = (1-rr/2)/(1+rr/2)*exp(nu2*rr^3)
        B = (1-rr/2)*(1+rr/2) + B2*rr^2
        
        Rgm, dtR = Rgmf(yn, X, Osb) #isoradial
        Xobi = X/Rgm #isoradial coordinate
        Xob = Xobi*B/enu #isotropic x to be referenced with rr
        #Xob = Xobi*enu/B #isotropic M/bar{R} to be referenced with rr 
        
        if rr > maxr
            maxr = rr
        end
        if rr < maxr
            if rr < Xob/0.95
                hit = false
                break
            end
        end
        
        if rr > Xob*0.95
            level = 2.0^6
            #level *= 2.0
        end

        #rr += hi
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
    #mult = 1.0
    
    t1 = (ep + Lz*w)
    #t1 = 1.0
    
    sq2 =  (ep + Lz*w)^2 - C*(enu^4)*(Xob^2)/(B^2)/(Rg^2)
    if sq2 < 0; sq2 = 0.0; end
    t2 = cosg*sqrt(sq2)
    #t2 = 0.0
    
    sq3 = C - Lz^2*csc(theta)^2
    if sq3 < 0; sq3 = 0.0; end
    #if sq3 < 0; sq3 = -sq3; end
    t3 = psign*sing*enu^2/B*(Xob/Rg)*sqrt(sq3)
    #t3 = 0.0
    
    cosa = 1 + mult*(-t1 + t2 + t3)
    
    if cosa < 0
        hit = false
    end
    cosa = clamp(cosa, 0.0, 1.0)
    
    #return rns, yns, zns, ers, lvs
    return time, phi, theta, Xob, hit, cosa
end

#println(bender3(0,6.5,sini,X,Osb,beta,quad,wp,Rg))
#println(bender3(0,8,sini,X,Osb,beta,quad,wp,Rg))

#IF photon path i.e. integrator debugger
if false

    #i=pi/2.01
    xpoint = 5.0
    ypoint = 0.004999999999999005
    #xpoint2 = 7.4374
    xpoint2 = 5.0
    ypoint2 = 0.20509999999999806

    #i=pi/4
    #xpoint = 5.0
    #ypoint = 0.05
    #xpoint2 = 0.5
    #ypoint2 = 6.5

    #i=0.05
    #xpoint = 3.0
    #ypoint = -3.0
    #xpoint2 = 3.0
    #ypoint2 = 3.0

    
rns, yns, zns, ers, lvs = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)

println(pi/2-yns[end])
p1 = plot(rns, (pi/2-yns), "b-")
rns, yns, zns, ers, lvs = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p1 = oplot(rns, (pi/2-yns), "r--")
println(pi/2-yns[end])
    
rns, yns, zns, ers, lvs = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p2 = plot(rns, zns, "b-")
rns, yns, zns, ers, lvs = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p2 = oplot(rns, zns, "r-")

rns, yns, zns, ers, lvs = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p3 = plot(rns, ers, "b--")
rns, yns, zns, ers, lvs = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p3 = oplot(rns, ers, "r--")

rns, yns, zns, ers, lvs = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p4 = plot(rns, lvs, "b-")
rns, yns, zns, ers, lvs = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p4 = oplot(rns, lvs, "r-")

t = Table(4,1)
t[1,1] = p1
t[2,1] = p2
t[3,1] = p3
t[4,1] = p4    
    #pval, tval, hit, cosa = bender3(3.0, 1.0, sini)

    #p = plot(rns, yns)
    

    
end


#IF compute image
if true
######################
######################
#grid setup
xmin = -10
xmax = -xmin +0.005
ymin = -10
ymax = -ymin +0.005

Ny = 101
Nx = 101
Nt = 1

xmid = round(Int,Nx/2)
ymid = round(Int,Ny/2)
hdata3 = zeros(Ny, Nx)

dx = (xmax-xmin)/(Nx-1)
dy = (ymax-ymin)/(Ny-1)
    
#x_grid = linspace(xmin, xmax, Nx)
#y_grid = linspace(ymin, ymax, Ny)

x_grid = collect(xmin:dx:xmax)
y_grid = collect(ymin:dy:ymax)   

Times = zeros(Ny,Nx)
Phis = zeros(Ny,Nx)
Thetas = zeros(Ny,Nx)
hits = zeros(Ny,Nx)
cosas = zeros(Ny,Nx)
Xs = zeros(Ny,Nx)

#Get middle point
time, phi, theta, Xob, hit, cosa = bender3(x_grid[xmid], y_grid[ymid], sini,
                                     X, Osb,
                                     beta, quad, wp, Rg)

Times[ymid, xmid] = time
Phis[ymid, xmid] = phi
Thetas[ymid, xmid] = theta


#cenval_t = theta
#cenval_p = phi
#midval_t = cenval_t
#midval_p = cenval_p

println("Computing image...")
tic()
#Upper region
for j = ymid:Ny
    #upper right corner
    if mod(j, 20) == 0; println("u: $j"); end
    #last_theta = midval_t
    #last_phi   = midval_p
    for i = xmid:Nx
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa
        #if hit
        #    last_theta = theta
        #    last_phi   = phi
        #end
        #if i == xmid
        #    midval_t = theta
        #    midval_p = phi
        #end
        if !hit
            break
        end
    end

    #upper left corner
    #last_theta = midval_t
    #last_phi   = midval_p
    for i = xmid:-1:1
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        #if hit
        #    last_theta = theta
        #    last_phi   = phi
        #end
        if !hit
            break
        end
    end
end

#Lower region
#midval_t = cenval_t
#midval_p = cenval_p
for j = ymid:-1:1
    if mod(j, 20) == 0; println("l: $j"); end
    #lower right corner
    #last_theta = midval_t
    #last_phi   = midval_p
    for i = xmid:Nx
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        #if hit
        #    last_theta = theta
        #    last_phi   = phi
        #end
        #if i == xmid
        #    midval_t = theta
        #    midval_p = phi
        #end
        if !hit
            break
        end
    end

    #lower left corner
    #last_theta = midval_t
    #last_phi = midval_p
    for i = xmid:-1:1
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                                   X, Osb,
                                                   beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        #if hit
        #    last_theta = theta
        #    last_phi   = phi
        #end
        if !hit
            break
        end
    end
end
toc()

end #if -tru


#IF interpolation
if true
    
p00 = plot2d(Times, x_grid, y_grid)
p11 = plot2d(Phis, x_grid, y_grid)
p21 = plot2d(Thetas, x_grid, y_grid)
p31 = plot2d(hits, x_grid, y_grid)

#interpolate into dense grid
println("interpolating into dense grid...")
#method = InterpQuadratic
method = InterpLinear
extrapolate = BCnearest

#dx = round(x_grid[2]-x_grid[1], 6)
#dy = round(y_grid[2]-y_grid[1], 6)
Xrange = xmin:dx:xmax
Yrange = ymin:dy:ymax

Times = Times .- Times[ymid, xmid]    
time_interp = CoordInterpGrid((Yrange, Xrange), Times,
                              extrapolate, InterpLinear)

phi_interp_sin = CoordInterpGrid((Yrange, Xrange), sin(Phis),
                             extrapolate, InterpLinear)
phi_interp_cos = CoordInterpGrid((Yrange, Xrange), cos(Phis),
                             extrapolate, InterpLinear)
        
#phi_interp = CoordInterpGrid((Yrange, Xrange), Phis,
#                             extrapolate, InterpLinear)
                             #extrapolate, InterpQuadratic)

theta_interp = CoordInterpGrid((Yrange, Xrange), Thetas,
                               extrapolate, method)

Xs_interp = CoordInterpGrid((Yrange, Xrange), Xs,
                             extrapolate, method)

cosa_interp = CoordInterpGrid((Yrange, Xrange), cosas,
                             extrapolate, method)

hits_interp = CoordInterpGrid((Yrange, Xrange), hits,
                             extrapolate, InterpLinear)



Ny_dense = 1000
Nx_dense = 1000
Times_dense = zeros(Ny_dense, Nx_dense)
Phis_dense = zeros(Ny_dense, Nx_dense)
Thetas_dense = zeros(Ny_dense, Nx_dense)
hits_dense = zeros(Ny_dense, Nx_dense)
    
img = zeros(Ny_dense, Nx_dense)
img2 = zeros(Ny_dense, Nx_dense) #debug array
img3 = zeros(Ny_dense, Nx_dense) #debug array
img4 = zeros(Ny_dense, Nx_dense) #debug array
        
Flux = zeros(Ny_dense, Nx_dense)
Reds = zeros(Ny_dense, Nx_dense)

painter = chess_board

x_grid_d = linspace(xmin, xmax, Nx_dense)
y_grid_d = linspace(ymin, ymax, Ny_dense)

    

function radiation(Ir,
                   x,y,
                   phi, theta, cosa,
                   X, Xob, Osb, sini)

    nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
    B2    = beta
    zeta2 = beta*(3*0.5*(3*cos(theta)^2-1)/4-1/3)
    Rgm, dR = Rgmf(theta, X, Osb)
    #Xob = X/Rgm

    enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
    B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
    ezeta = (1-Xob/2)*(1+Xob/2) + zeta2*Xob^2

    #println("xobi:",Xob*enu/B)
    
    C = (x^2 + y^2)
    Lz = x*sini
    w = wp*Xob^3*(1-3*Xob)

    fa = (B/enu/ezeta)*dR/Rgm
    #fa = dR/Rgm
    
    cosg = 1/sqrt(1 + fa^2)
    sing = fa*cosg

    #if w == 0
    if false
    ####
    #vphi = (Rgm/enu)*(2*pi*fs - w)*sin(theta)
    #vphi = Rgm*B*(1/enu^2)*sin(theta)*(2pi*fs - w) #isotropic zamo WRONG
    vphi = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    b = R*vphi/c

    #vw = Rgm*B*(1/enu^2)*sin(theta)*w #isotropic space vel WRONG
    vw = Rgm*(1/enu)*sin(theta)*w #isoradial space vel
    bp = R*vw/c
        
    gamma = 1/sqrt(1 - b^2)
    cosi = sqrt(1-sini^2)
    sina = sqrt(1-cosa^2)
    cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
        
    cosz = -sina*sini*sin(phi)/sqrt(1-cospsi^2)
 
    delta = (1/gamma)/(1 - b*cosz)
    EEd = delta*enu*(1 + cosz*bp)

    ###
    else
    
    #vz = Rgm*B*(1/enu^2)*sin(theta)*(2pi*fs - w) #isotropic zamo
    vz = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    bz = R*vz/c

    gamma = 1/sqrt(1-bz^2)
    delta = 1/gamma/(1 + Lz*(2pi*fs)/(G*M/c^2))
    EEd = delta*enu
    end

    
    dS = (Rgm)^2*sin(theta)*sqrt(1 + fa^2)
    cosap = cosa *delta
    dOmega = dS*cosap
    
    dF = (EEd^3)*dOmega*Ir(cosap)
    
    return dF, EEd
    #return gamma, gamma2
    #return dS, dOmega
end

    dxx = dx
    dyy = dy

    dxx = 1.0*(x_grid_d[2] - x_grid_d[1])
    dyy = 1.*(y_grid_d[2] - y_grid_d[1])
    dxdy = dxx*dyy*X^2

phi_interp_atan(y,x) = atan2(phi_interp_sin[y,x],phi_interp_cos[y,x])

##############################
for j = 1:Ny_dense
    y = y_grid_d[j]
    for i = 1:Nx_dense
        
        # phi & theta
        
        
        x = x_grid_d[i]
        phi = phi_interp_atan(y,x)
        #phi = phi_interp[y,x]
        theta = theta_interp[y,x]
        Xob = Xs_interp[y,x]
        time = time_interp[y,x]
        
        #test if we hit the surface
        hit = hits_interp[y,x]
        hiti = round(Int,hit - 0.1)

        if hiti > 0
            #solid angle
            ####
            
            #squares
            function polyarea(x,y,dxx,dyy,phi0,the0)
                              
                #image plane
                x1 = x - dxx/2
                x2 = x + dxx/2
                y1 = y - dyy/2
                y2 = y + dyy/2
            
                #surface
                phi1 = phi_interp_atan(y1, x1)
                phi2 = phi_interp_atan(y1, x2)
                phi3 = phi_interp_atan(y2, x2)
                phi4 = phi_interp_atan(y2, x1)
                the1 = theta_interp[y1, x1]
                the2 = theta_interp[y1, x2]
                the3 = theta_interp[y2, x2]
                the4 = theta_interp[y2, x1]
                
                #
                #phi = (mod2pi(phi1)+mod2pi(phi2)+mod2pi(phi3)+mod2pi(phi4))/4
                #theta = (mod2pi(the1)+mod2pi(the2)+mod2pi(the3)+mod2pi(the4))/4
                #parea = area_sphere([phi1, phi2, phi3, phi4],
                #                    [the1, the2, the3, the4],
                #                    Rq, ecc)
                parea = area_sphere(phi0,the0,
                                    [phi1, phi2, phi3, phi4],
                                    [the1, the2, the3, the4],
                                    Rq, ecc)


                
                return parea
            end
            #triangles
            function polyarea2(x,y,dxx,dyy)
                              
                #image plane
                x1 = x - dxx
                x2 = x + dxx
                y1 = y - dyy/2
                y2 = y + dyy/2
            
                #surface
                phi1 = phi_interp_atan(y1, x1)
                phi2 = phi_interp_atan(y1, x2)
                phi3 = phi_interp_atan(y2, x)
                the1 = theta_interp[y1, x1]
                the2 = theta_interp[y1, x2]
                the3 = theta_interp[y2, x]

                #phi = (mod2pi(phi1)+mod2pi(phi2)+mod2pi(phi3)+mod2pi(phi4))/4
                #theta = (mod2pi(the1)+mod2pi(the2)+mod2pi(the3)+mod2pi(the4))/4
                #parea = area_sphere([phi1, phi2, phi3],
                #                    [the1, the2, the3],
                #                    Rq, ecc)
                parea = area_sphere_tri([phi1, phi2, phi3],
                                        [the1, the2, the3],
                                        Rq, ecc)
                
                return parea
            end

            
            #center value
            Sig = dxdy/polyarea(x,y,dxx,dyy,phi,theta)
            #Sig = polyarea(x,y,dxx,dyy)
            #corners
            #Sig1 = dxdy/polyarea(x-dxx,y-dyy,dxx,dyy)
            #Sig2 = dxdy/polyarea(x+dxx,y-dyy,dxx,dyy)
            #Sig3 = dxdy/polyarea(x-dxx,y+dyy,dxx,dyy)
            #Sig4 = dxdy/polyarea(x+dxx,y+dyy,dxx,dyy)
            #Sig = (Sig+Sig1+Sig2+Sig3+Sig4)/5.0
            
            
        Phis_dense[j,i] = phi
        Thetas_dense[j,i] = theta
        Times_dense[j,i] = time
            
        #chess board
        img[j,i] = painter(phi, theta)
            
        #mu = sqrt(1-sini^2)*cos(theta) + sini*sin(theta)*cos(phi)
        #img2[j,i] = mu

        #radiation
        cosa = cosa_interp[y,x]
        #if 0 < cosa < 1

        nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
        B2    = beta
        zeta2 = beta*(3*0.5*(3*cos(theta)^2-1)/4-1/3)
        
        Rgm, dR = Rgmf(theta, X, Osb)
        #cosa = cosalpha(x, y, sini, Rgm, dR,
        #                X/Rgm, nu2, B2, zeta2, wp, theta, Rg)

                
        img2[j,i] = cosa


        ##################################
        ##################################
        #approximative cosa test
        #if false
        if true        

            nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
            B2    = beta
            zeta2 = beta*(3*0.5*(3*cos(theta)^2-1)/4-1/3)
            Rgm, dR = Rgmf(theta, X, Osb)
            #Xob = X/Rgm
            
            enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
            B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
            ezeta = (1-Xob/2)*(1+Xob/2) + zeta2*Xob^2

            w = wp*Xob^3*(1-3*Xob)*G*M/c^2
            
            vphi = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
            b = R*vphi/c

            #vw = Rgm*B*(1/enu^2)*sin(theta)*w #isotropic space vel WRONG
            vw = Rgm*(1/enu)*sin(theta)*w #isoradial space vel
            bp = R*vw/c
        
            gamma = 1/sqrt(1 - b^2)
            cosi = sqrt(1-sini^2)
            sina = sqrt(1-cosa^2)
            cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
            cosz = -sina*sini*sin(phi)/sqrt(1-cospsi^2)

            ut = U/Rgm
            #xt = X/Rgm

            sqrta(xxx) = sqrt(abs(xxx))
            #cospsi = sqrt(1-sini^2)*cos(theta)+sini*sin(theta)*cos(phi)
            #cosaa = 1- (1-cospsi)*(1-Xob/2)/(1+Xob/2)*exp(nu2*(Xob)^3)
            #cosaa = 1 - (1-cospsi)*(1-ut)

            img4[j,i] = 1 - (1-cospsi)*(enu^2)#non-rotating reference approx
            
            cosaa = 1 - (1-cospsi)*(enu^2)/(1+cosz*bp)^3 #*(1-bp^2)^2
            #cosaa = 1 - (1-cospsi)*(enu^2)/(1+cosz*bp)^2.0 /((1-bp^2)^2)
            #cosaa = 1 - (1-cospsi)*(enu^2)*(1 + cosz*bp)^2 #Sul

            #cosaa = 1 - (1-cospsi)*(enu^2)/(1+cosz*(b+bp))^2 #/((1-bp^2)/(1+bp*b))^2
            
            #initially ingoing photons
            #bimpact = Rgm*sqrta(1-cosaa^2)/sqrta(1-ut)
            #bimpact = Rgm*sqrta(1-cosaa^2)/( (1-xt/2)/(1+xt/2) )
            #if bimpact < Rgm/sqrta(1-ut) && acos(cosaa) > pi/2
                
                #    rcsqrt = -(-9*bimpact^2*Rgm*ut + sqrta(12*bimpact^6 - 81*bimpact^4*Rgm^2*ut^2))
                #    rc = -(-2*3^(1/3)*bimpact^2 + 2^(1/3)*rcsqrt^(2/3))/(6^(2/3)*rcsqrt^(1/3))
                #    dspi = sqrta(2/rc)*sqrta(Rgm-rc)/sqrta(1-ut*(Rgm/rc))
                #    ocosaa = cosaa
                #    cosaa = cos(dspi + acos(ocosaa))
                #cosaa = -10.0
            #end

        #Change from spherical system to oblate according to Morsink et al 2007
        fa = (1/sqrt(1-ut))*dR/Rgm
        #cf = dR/Rgm
        cosg = 1/sqrt(1 + fa^2)
        sing = fa*cosg

        #cosg=0.5
        #sing=0.5
            
        #sing = sqrt(1-cosg^2)
        cosd = (sqrt(1-sini^2)-cos(theta)*cospsi)/(sin(theta)*sqrt(1-cospsi^2))
        if abs(sin(phi)) < 1.0e-4
            cosb = cosaa*cosg
        else
            cosb = cosaa*cosg + sqrta(1-cosaa^2)*sing*cosd
        end
        #cosb = clamp(cosb, 0.,1.)

        #img3[j,i] = cosb
        #img3[j,i] = fa    
            #img3[j,i] = cosaa
            #img3[j,i] = cosz
            #img3[j,i] = cosz*bp
        end #if false/true for cosa
            
        #img3[j,i] = time    
        img3[j,i] = Sig

        # Solid angle
        
            
        ##################################
        ##################################

        #Radiation
        Ir(cosa) = 1.0 #isotropic beaming
            
        dF, dE = radiation(Ir,
                           x,y,
                           phi, theta, cosa,
                           X, Xob, Osb, sini)

        #if 0.79 < dE < 0.81
        #if dE < 0.79 || dE > 0.81
        Flux[j,i] = dF
        Reds[j,i] = dE
        #end
            
        end#hiti
        #end#cosa
    end
end

#Phis = Phis_dense
        #Thetas = Thetas_dense
        #println("mean:",mean(img3))
#img3 = img3 ./ mean(img3)
    
p0 = plot2d(Times_dense, x_grid_d, y_grid_d, 0,0,0,"Blues")
p1 = plot2d(Phis_dense, x_grid_d, y_grid_d)
p2 = plot2d(Thetas_dense, x_grid_d, y_grid_d)
p3 = plot2d(img, x_grid_d, y_grid_d)

p4 = plot2d(img2, x_grid_d, y_grid_d, 0, 0, 0, "Blues")
p5 = plot2d(img3, x_grid_d, y_grid_d, 0, 0.0, 1.2, "Blues")

p6 = plot(y_grid_d, img2[:,511],"k-", yrange=[-0.1, 1.1])
p6 = oplot(y_grid_d, img3[:,511], "r--")
p6 = oplot(x_grid_d, img2[501,:], "b-")
p6 = oplot(x_grid_d, img3[501,:], "r--")

rel_err(x1,x2) = (x1 .- x2) ./ x1
xslice = 501
yslice = 501

p6e = plot(y_grid_d, rel_err(img2[:,xslice],img3[:,xslice]), "k-",yrange=[-0.05, 0.05])
p6e = oplot(y_grid_d, rel_err(img2[yslice,:],img3[yslice,:]), "b-")
p6e = oplot(y_grid_d, rel_err(img2[yslice,:],img4[yslice,:]), "g", linestyle="dotted")

p6e = oplot(y_grid_d, zeros(length(y_grid_d)), "k",linestyle="dotted")
p6e = oplot(y_grid_d, img2[:,xslice]*0.02, "k",linestyle="solid")
p6e = oplot(y_grid_d, img2[yslice,:]*0.02, "k",linestyle="solid")
p6e = oplot(y_grid_d, img3[:,xslice]*0.02, "r",linestyle="dashed")
p6e = oplot(y_grid_d, img3[yslice,:]*0.02, "m",linestyle="dashed")

p6e2 = plot(img2[:,xslice], rel_err(img2[:,xslice],img3[:,xslice]), "k-",yrange=[-0.6, 0.3])
p6e2 = oplot(img2[:,xslice], rel_err(img2[yslice,:],img3[yslice,:]), "b-")
p6e2 = oplot(img2[:,xslice], rel_err(img2[yslice,:],img4[yslice,:]), "g",linestyle="dotted")
p6e2 = oplot(img2[:,xslice], zeros(length(y_grid_d)), "k",linestyle="dotted")

#p6e2 = plot2d(rel_err(img2,img3), x_grid_d, y_grid_d, 0,0,0,"Blues")

#p7 = plot2d(Flux, x_grid_d, y_grid_d, 0,0,0, "RdBu")
#p8 = plot2d(Reds, x_grid_d, y_grid_d, 0,0,0, "RdBu")
p7 = plot2d(Flux, x_grid_d, y_grid_d, 0,0,0, "Blues")
p8 = plot2d(Reds, x_grid_d, y_grid_d, 0,0,1, "Blues")




#####
# line profile
function line_prof(Flux, Reds)
    println("Computing line profile...")
    #for energ = linspace(1-0.01, 1+0.01, 40)

    xarr = Float64[]
    yarr = Float64[]

    Ny_dense, Nx_dense = size(Flux)
    
    energ = 1.0
    for jj = 1:Ny_dense, ii = 1:Nx_dense
        fy =  Flux[jj, ii]
        xii = Reds[jj, ii]
        
        if xii > 0
            push!(xarr, xii*energ)
            push!(yarr, fy)
        end
    end
    #end
    
    xind = sortperm(xarr)
    xarrs = xarr[xind]
    yarrs = yarr[xind]
    NN = length(xarrs)
    
    emin = minimum(xarrs)
    emax = maximum(xarrs)
    println("emin=$emin emax=$emax")
    Nr = 150
    es = collect(linspace(emin, emax, Nr))
    yy2 = zeros(Nr)
    
    xst = 1
    for ii = 2:Nr
        for jj = xst:NN
            if es[ii-1] <= xarrs[jj] < es[ii]
                yy2[ii] += yarrs[jj]
            elseif xarrs[jj] >= es[ii]
                xst = jj
                break
            end
        end
    end
    
    yy2 = yy2./maximum(yy2)
    #add start and end points to make smooth figure
    unshift!(es, es[1])
    unshift!(yy2, 0.0)
    push!(es, es[end])
    push!(yy2, 0.0)

    return es, yy2
end

es, yy2 = line_prof(Flux, Reds)
p9 = plot(es, yy2, "k-")
          #xlabel="E/E_0",
          #ylabel="Flux (arb)")
          #xrange=[0.7, 0.9])

end



#collect into table
tY = 3
tX = 4
table = Table(tY,tX)
images = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9]
imc = 1
for tty = 1:tY
    for ttx = 1:tX
        table[tty,ttx] = images[imc]
        imc += 1

        if imc > length(images)
            break
        end
    end
end
display(table)


#locate star edges for integration limits
println("Computing waveform...")
tic()

x1s = zeros(Int, Ny_dense)
x2s = zeros(Int, Ny_dense)
y1s = 0
y2s = 0

top = false
for j = 1:Ny_dense
    y = y_grid_d[j]

    left = false
    for i = 1:Nx_dense
        x = x_grid_d[i]
        
        hit = hits_interp[y,x]
        hiti = round(Int,hit - 0.1)

        if hiti == 1
            if left == false
                x1s[j] = i
                left = true
            else
                x2s[j] = i
            end

            if top == false
                y1s = j
                top = true
            else
                y2s = j
            end
        end
    end
end

println()
println("size (isotropic)")
println("xmin=",x_grid_d[minimum(x1s[y1s:y2s])]," xmax=",x_grid_d[maximum(x2s[y1s:y2s])])
println("ymin=",y_grid_d[y1s]," ymax=",y_grid_d[y2s])


Xob = Xs_interp[0,0]
#enu = (1-Xob/2)/(1+Xob/2)
#B=(1-Xob/2)*(1+Xob/2)
nu2   = beta/3.0 - quad*0.5*(-1.0)
B2    = beta
enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
corr = enu/B

println()
println("size (isoradial)")
println("xmin=",x_grid_d[minimum(x1s[y1s:y2s])]*corr," xmax=",x_grid_d[maximum(x2s[y1s:y2s])]*corr)
println("ymin=",y_grid_d[y1s]*corr," ymax=",y_grid_d[y2s]*corr)


#IF Spot animation
if false

    
########
function spot(t, phi, theta;
              stheta = pi/3, #spot latitude
              delta = 40*pi/180 #spot half-opening angle
              )

    #Vincenty's formula
    d = great_circle_dist(0.0, phi, stheta, theta)
    
    if abs(d) < delta
        return true
    end
    
    return false
end

img4 = zeros(Ny_dense, Nx_dense) #debug array
Ir(cosa) = 1.0 #isotropic beaming
Nt = 60
times = linspace(0, 1/fs, Nt)


for k = 1:Nt
    img4[:,:] = 0.0
    
    t = times[k]
    for j = y1s:y2s
        y = y_grid_d[j]
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            theta = theta_interp[y,x]

            #rotate star
            phi = phi_interp_atan(y,x)
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)
      

            img4[j,i] = painter(phi, theta)/10

            inside = spot(0.0, phi, theta)
            
            if inside
            #if true    
                time = time_interp[y,x]
                Xob = Xs_interp[y,x] 
                cosa = cosa_interp[y,x]
                dF, dE = radiation(Ir,
                                   x,y,
                                   phi, theta, cosa,
                                   X, Xob, Osb, sini)

                
                #img4[j,i] = painter(phi, theta)
                img4[j,i] += dF
                
            end #if inside            
        end# for x
    end#for y

    p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 0, "Blues")
    #p10 = plot2d(img4, x_grid_d, y_grid_d, 0, -1., 1., "RdBu")
    display(p10)
    
end#for t


toc()
p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 0, "Blues")


end #if false true -waveform
