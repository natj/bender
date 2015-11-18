#Area of a polygon on a sphere/spheroid
function area_sphere(lons, cols, Rarea=1, ecc=1)

    N = length(cols)
    
    #transform to authalic sphere if we have spheroid
    if ecc < 1.0
        for i = 1:N
            cols[i] = authalic_lat(cols[i], ecc, qp)
        end
    end

    #sigs = zeros(N)
    #for i = 1:N-1
    #    sigs[i] = great_circle_dist(lons[i], lons[i+1], cols[i], cols[i+1]) 
    #end
    #sigs[N] = great_circle_dist(lons[N], lons[1], cols[N], cols[1]) 
    #area = Rarea^2 * (sum(abs(sigs)) - (N-2)*pi)
    #return area
    
    #area = (lons[2] - lons[N])*cos(cols[1])
    #for i = 2:N-1
    #    area += (lons[i+1] - lons[i-1])*cos(cols[i])
    #end
    #area += (lons[1] - lons[N-1])*cos(cols[N])
    #area = (area*R^2)/2
    #return area
    
    
    function azimuth(lon1, lon2, col1, col2)
        lat1 = col1 - pi/2
        lat2 = col2 - pi/2
        
        f1 = cos(lat2) * sin(lon2-lon1)
        f2 = cos(lat1) * sin(lat2)
        f3 = sin(lat1) * cos(lat2) * cos(lon2-lon1)
        az = atan2(f1, f2-f3)
        return mod2pi(az)
    end

    #lats = cols - pi/2
    #lat0 = 0
    col0 = pi/2
    lon0 = 0
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

lats = [-90:30:90,60:-30:-60] .- 90
lons = [zeros(7),90*ones(5)]
a = area_sphere(deg2rad(lons), deg2rad(lats))/(4*pi)
println(a)



lats2=[-90:1:90,89:-1:-89] .- 90
lons2=[zeros(181),90*ones(179)]
a = area_sphere(deg2rad(lons2), deg2rad(lats2))/(4*pi)
println(a)
