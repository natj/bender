######################
######################

#Spherical radial function
function Rgmf2(theta, X, Osb)
    return 1.0, 0.0
end

#simplified radial function
function Rgmf4(theta, X, Osb)
    Rgm = (1-0.15*cos(theta)^2)
    dtR = 0.3*sin(theta)*cos(theta)

    return Rgm, dtR
end


#Radial function by AlGendy & Morsink 2014
#Circumferential (isoradial) radius
function Rgmf(theta, X, Osb)
    const o20 = -0.788
    const o21 = 1.030

    #Radial function
    req=1.0
    o2=(o20+o21*X)*Osb^2
    Rgm = req*(1+o2*cos(theta)^2)

    #println("Osb = $Osb")
    #println("o2 = $o2")
    #println("dr = ", -2*Osb^2*(o20+o21*X))

    #derivative dR/dtheta
    dtR = -2*Osb^2*(o20+o21*X)*cos(theta)*sin(theta) #from AlGendy & Morsink 2014

    return Rgm, dtR
end


######################
######################
#Radial function by Cadeau 2007
#Circumferential radius
function Rgmf3(theta, X, Osb)
    
    angvel = 2*pi*fs
    zeta = G*M/(R*c^2)
    eps = (angvel^2 * R^2)/(zeta*c^2)
    #println("zeta: $zeta  eps: $eps")
    #eps2 = (angvel^2 * R^3)/(G*M)
    #eps = 0.339
    
    a0 = -0.18*eps + 0.23*zeta*eps - 0.05*eps^2
    a2 = -0.39*eps + 0.29*zeta*eps + 0.13*eps^2
    a4 =  0.04*eps - 0.15*zeta*eps + 0.07*eps^2

    #CFL Quark Star
    #    a0 = -0.26*eps+0.50*zeta*eps-0.04*eps*eps 
    #    a2 = -0.53*eps+0.85*zeta*eps+0.06*eps*eps
    #    a4 = 0.02*eps-0.14*zeta*eps+0.09*eps*eps
    #end
    
    #Legendre polynomials
    p0 = 1.0
    p2 = 0.25*(1.0 + 3.0*cos(2.0*theta))
    p4 = 0.015625*(9.0 + 20.0*cos(2.0*theta) + 35.0*cos(4.0*theta))
    
    sradius = (1.0 + a0*p0 + a2*p2 + a4*p4)

    #derivative
    dp2 = -3.0*sin(theta)*cos(theta)
    dp4 = -5/16*(2*sin(2*theta) + 7*sin(4*theta))
    dtR = a2*dp2 + a4*dp4

    #return 1.0, 0.0
    return sradius, dtR
end



###################
# Trigonometric functions


#Eccentricity for area projections
eRe, dR = Rgmf(pi/2, X, Osb) #equatorial radius
eRp, dR = Rgmf(0, X, Osb) #pole radius

const fe = (eRe - eRp)/eRe
const ecc = sqrt(1 - (eRp/eRe)^2) #eccentricity
if ecc != 0.0
    const qp = 1 + (1-ecc^2)/(2*ecc)*log((1+ecc)/(1-ecc)) #authalic pole latitude
    const Rq = eRe*sqrt(qp/2) #authalic radius
else
    const qp = 0.0
    const Rq = 1.0
end
println("ecc = $ecc Rq = $Rq fe = $fe")


#authalic colatitude 
function authalic_lat(colat, ecc, qp)

    lat = pi/2 - colat
    sinl = sin(lat)
    q = (1-ecc^2)*(sinl/(1 - (ecc*sinl)^2) - log((1 - ecc*sinl)/(1 + ecc*sinl))/2/ecc)
    alat = asin(q/qp)
    acolat = pi/2 - alat
    
    return acolat
end

#Great circle distance (on a sphere)
function great_circle_dist(lon1, lon2, col1, col2)

    lat1 = pi/2 - col1
    lat2 = pi/2 - col2
    #lat1 = col1
    #lat2 = col2
    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)

    #law of cosines
    #return acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon))
    
    #haversine formula
    #return 2.0*asin(sqrt(sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2))

    #vincenty's
    xx = (cos(lat2)*sin(dlon))^2 + (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))^2
    yy = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
    return atan2(sqrt(xx), yy)
end

function great_circle_dist_rel(lon1, lon2, col1, col2, gamma)

    #lat1 = pi/2 - col1
    #lat2 = pi/2 - col2
    lat1 = col1
    lat2 = col2
    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)

    #relativistic law of cosines 
    gamma1 = gamma
    gamma2 = gamma
    norm1=cos(lat1)^2 + (sin(lat1)^2)*(cos(lon1)^2 + (gamma1*sin(lon1))^2)
    norm2=cos(lat2)^2 + (sin(lat2)^2)*(cos(lon2)^2 + (gamma2*sin(lon2))^2)
    t1=cos(lat1)*cos(lat2) + sin(lat1)*sin(lat2)*(cos(lon1)*cos(lon2) + gamma1*gamma2*sin(lon1)*sin(lon2))
    return acos(t1/sqrt(norm1)/sqrt(norm2))


    #relativistic haversine
    gamma1 = gamma
    gamma2 = gamma
    t1=(-cos(lat2)*cos(lon1)*sin(lat1) + cos(lat1)*cos(lon2)*sin(lat2))^2

    t2=(gamma1*cos(lat2)*sin(lat1)*sin(lon1) - gamma2*cos(lat1)*sin(lat2)*sin(lon2))^2
    
    t3 = (sin(lat1)*sin(lat2)*(-gamma1*cos(lon2)*sin(lon1) + gamma2*cos(lon1)*sin(lon2)))^2

    #t3=(gamma1*cos(lon2)*sin(lat1)*sin(lat2)*sin(lon1) - 
    #    gamma2*cos(lon1)*sin(lat1)*sin(lat2)*sin(lon2))^2

    norm1=cos(lat1)^2 + (sin(lat1)^2)*(cos(lon1)^2 + (gamma1*sin(lon1))^2)
    norm2=cos(lat2)^2 + (sin(lat2)^2)*(cos(lon2)^2 + (gamma2*sin(lon2))^2)
    return asin(sqrt(t1+t2+t3)/sqrt(norm1)/sqrt(norm2))

    
    #haversine formula
    #return 2.0*asin(sqrt(sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2))

    #vincenty's
    xx = (cos(lat2)*sin(dlon))^2 + (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon))^2
    yy = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
    #return atan2(sqrt(xx), yy)
end
#Great circle distance (on an ellipsoid)
#uses general Vincenty's formula for oblate spheroids
function great_circle_dist2(lon1, lon2, col1, col2)

    lat1 = pi/2 - col1
    lat2 = pi/2 - col2
    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)

    #Vincencty's formula on an ellipsoid (wiki)
    u1 = atan((1-fe)*tan(lat1))
    u2 = atan((1-fe)*tan(lat2))
    lambda = dlon
    err = 1.0

    omega = 1.
    cospalpha = 1.
    sinomega = 1.
    cosomega = 1.
    costomega = 1.
    cospomega = 1.
        
    while err > 1e-8
        sinomega = sqrt( (cos(u2)*sin(lambda))^2 + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(lambda))^2)
        cosomega = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(lambda)
        omega = atan2(sinomega, cosomega)

        sinalpha = cos(u1)*cos(u2)*sin(lambda)/sinomega
        cospalpha = 1 - sinalpha^2
        costomega = cosomega - 2*sin(u1)*sin(u2)/cospalpha

        Cf = (fe/16)*cospalpha*(4+fe*(4 - 3*cospalpha))
        lambdan = dlon + (1-Cf)*fe*sinalpha*(omega + Cf*sinomega*(costomega) + Cf*cosomega*(-1 + 2*cospomega))

        err = (lambdan-lambda)/lambdan
        lambda = lambdan
    end

    us2 = cospalpha*(eRe^2 - eRp^2)/eRp^2
    A = 1 + us2/16384*(4096 + us2*(-768 + us2*(320 - 175*us2)))
    B = us2/1024*(256 + us2*(-128 + us2*(74 - 47*us2)))
    domega = B*sinomega*(costomega + 0.25*B*(cosomega*(-1 + 2*costomega^2)-(1/6)*B*costomega*(-3 + 4*sinomega^2)*(-3 + 4*costomega^2)))

    s = eRp*A*(omega - domega)

    return s
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
    
    function azimuth(lon1, lon2, col1, col2)
        lat1 = col1 - pi/2
        lat2 = col2 - pi/2
        
        f1 = cos(lat2) * sin(lon2-lon1)
        f2 = cos(lat1) * sin(lat2)
        f3 = sin(lat1) * cos(lat2) * cos(lon2-lon1)
        az = atan2(f1, f2-f3)
        return mod2pi(az)
    end

    col0 = pi/2 - the0 
    lon0 = phi0

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

#Area on top of a sphere using lambert projection to cartesian space
function area_sphere_lambert(phi0, theta0, phis, thetas, Rarea=1, ecc=0.0)
    N = length(thetas)

    #for i = 1:N
    #    phis[i] = mod2pi(phis[i])
    #end
    
    #transform to authalic sphere if we have spheroid
    if ecc != 0.0
        for i = 1:N
            thetas[i] = authalic_lat(thetas[i], ecc, qp)
        end
        theta0 = authalic_lat(theta0, ecc, qp)
    end

    lphis = pi/2 - thetas
    lphi0 = pi/2 - theta0
    llambdas = phis
    llambda0 = phi0

    kprime = sqrt(2 ./ (1 .+ sin(lphi0) .* sin(lphis) .+ cos(lphi0) .* cos(lphis) .* cos(llambdas-llambda0)))
    xs = kprime .* cos(lphis) .* sin(llambdas-llambda0)
    ys = kprime .* (cos(lphi0) .* sin(lphis) .- sin(lphi0) .* cos(lphis) .* cos(llambdas-llambda0))

    xys = collect(zip(xs, ys))
    xys2 = append!(xys[2:end], [xys[1]])
    polygon = collect(zip(xys, xys2))
    area = 0.5 * abs(sum(
        [x0*y1 - x1*y0 for ((x0,y0), (x1,y1)) in polygon]))

    return area*4*pi*Rarea
end

#Area on top of a sphere using lambert projection to cartesian space
#using sin & cos of phi for numerical stability
function area_sphere_lambert2(phi0, theta0, sinphis, cosphis, thetas, Rarea=1, ecc=0.0)
    N = length(thetas)

    #for i = 1:N
    #    phis[i] = mod2pi(phis[i])
    #end

    #println()
    #println("phi0:$phi0 the0:$theta0 $sinphis $cosphis $thetas")
    
    #transform to authalic sphere if we have spheroid
    if ecc != 0.0
        for i = 1:N
            thetas[i] = authalic_lat(thetas[i], ecc, qp)
        end
        theta0 = authalic_lat(theta0, ecc, qp)
    end

    lphis = pi/2 - thetas
    lphi0 = pi/2 - theta0
    #llambdas = phis
    llambda0 = phi0

    #kprime = sqrt(2 ./ (1 .+ sin(lphi0) .* sin(lphis) .+ cos(lphi0) .* cos(lphis) .* cos(llambdas-llambda0)))
    #xs = kprime .* cos(lphis) .* sin(llambdas-llambda0)
    #ys = kprime .* (cos(lphi0) .* sin(lphis) .- sin(lphi0) .* cos(lphis) .* cos(llambdas-llambda0))

    coslml0 = sin(llambda0)*sinphis .+ cos(llambda0)*cosphis
    sinlml0 = cos(llambda0)*sinphis .- sin(llambda0)*cosphis
    
    #sqkprime = 2 ./ (1 .+ sin(lphi0) .* sin(lphis) .+ cos(lphi0) .* cos(lphis) .* coslml0)
    #println(sqkprime)

    #println("lphi0: $lphi0")
    #println("lphi: $lphis")
    #println("acos: $(coslml0)")
    
    kprime = sqrt(2 ./ (1 .+ sin(lphi0) .* sin(lphis) .+ cos(lphi0) .* cos(lphis) .* coslml0))
    xs = kprime .* cos(lphis) .* sinlml0
    ys = kprime .* (cos(lphi0) .* sin(lphis) .- sin(lphi0) .* cos(lphis) .* coslml0)

    
    xys = collect(zip(xs, ys))
    xys2 = append!(xys[2:end], [xys[1]])
    polygon = collect(zip(xys, xys2))
    area = 0.5 * abs(sum(
        [x0*y1 - x1*y0 for ((x0,y0), (x1,y1)) in polygon]))

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

