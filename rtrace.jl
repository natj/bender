#Compute image using polar coordinate system
include("bender.jl")


#grid setup
Nrad = 100
Nchi = 100

rmin = 0.0
rmax = 12.0

dchi_edge = 0.01
chimin = 0.0 - dchi_edge
chimax = 2.0pi + dchi_edge
#chimin = -pi
#chimax = pi


#get radius limits for the image
#chis = [0.0, pi/2, pi, 3pi/2]
Nedge = 15
chis = linspace(0.0, 1.0, Nedge)*2pi
rlims = zeros(length(chis))
for i = 1:length(chis)
    chii = chis[i]
    
    rmini = rmin
    rmaxi = rmax
    rmid = 0.0
    
    Nbi = 10
    N = 1
    while N <= Nbi
        rmid = (rmini + rmaxi)/2.0

        time, phi, theta, Xob, hit, cosa = bender3p(rmid, chii, sini,
                                                    X, Osb, beta, quad, wp, Rg)
        
        if hit
            rmini = rmid
        else
            rmaxi = rmid
        end
        N += 1
    end
    rlims[i] = rmid
end

#set new maximum rad limit
#rscale = maximum(rlims)
rmax = maximum(rlims)*1.02
println("max edge: $rmax")

#Create edge function to get the exact shape of the outline
#method = Gridded(Linear())
#edge_interp    = interpolate((chis), rlims, method)
edge_interp_raw = interpolate(rlims, BSpline(Quadratic(Periodic())), OnCell())
edge_interp(x) = edge_interp_raw[(Nedge-1)*x/2pi + 1.0]

#p = plot(chis/2pi, rlims, "rx")
#p = oplot(linspace(0.0, 1.0, 100), Float64[edge_interp(2pi*x) for x in linspace(0.0, 1.0, 100)], "b-")
#display(p)

#rad & chi grids
#equidistant
#chi_grid = collect(linspace(chimin, chimax, Nchi))
#rad_grid = collect(linspace(rmin, rmax, Nrad))

#equidistant
#chi_grid = collect(linspace(chimin, chimax, Nchi-2))

#weighted
chi_diffs = 0.8 + sin(collect(linspace(0.0, 2pi, Nchi-3))).^2
unshift!(chi_diffs, 0.0)
chi_grid = chimin .+ (chimax-chimin)*cumsum(chi_diffs)/sum(chi_diffs)

unshift!(chi_grid, chi_grid[1] - dchi_edge)
push!(chi_grid, chi_grid[end] + dchi_edge)


rad_diffs = 1 ./ exp(linspace(0.0, 2.0, Nrad-1).^2)
rad_grid = rmax * cumsum(rad_diffs) / sum(rad_diffs)
unshift!(rad_grid, 0.0)

#differential elements
drad = diff(rad_grid)
dchi = diff(chi_grid)



Times = zeros(Nrad, Nchi)
Phis = zeros(Nrad, Nchi)
Thetas = zeros(Nrad, Nchi)
hits = zeros(Nrad, Nchi)
cosas = zeros(Nrad, Nchi)
Xs = zeros(Nrad, Nchi)

#edge = zeros(Nchi)

println()
println("Computing image...")
tic()
for i = 1:Nchi
    chi = chi_grid[i]
    
    if mod(i,10) == 0; println("chi: $(round(chi/2pi,2))"); end
    
    for j = 1:Nrad
        rad = rad_grid[j]

        #if iseven(i) && rad < 0.66*rmax
        #    Times[j,i] = Times[j,i-1]
        #    Phis[j,i] = Phis[j,i-1]
        #    Thetas[j,i] = Thetas[j,i-1]
        #    Xs[j,i] = Xs[j,i-1]
        #    hits[j,i] = hits[j,i-1]
        #    cosas[j,i] = cosas[j,i-1]
        #    hit = convert(Bool, hits[j,i])
        #else
        time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                    X, Osb,
                                                    beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa
        #end

        #edge[i] = rad
        if !hit; break; end
    end
end
toc()

print("interpolating into dense grid...")
method = Gridded(Linear())
#method = Gridded(Constant())

time0 = Times[1,1]
Times = Times .- time0

time_interp    = interpolate((rad_grid, chi_grid), Times     ,method)
phi_interp_sin = interpolate((rad_grid, chi_grid), sin(Phis) ,method)
phi_interp_cos = interpolate((rad_grid, chi_grid), cos(Phis) ,method)
theta_interp   = interpolate((rad_grid, chi_grid), Thetas    ,method)
Xs_interp      = interpolate((rad_grid, chi_grid), Xs        ,method)
cosa_interp    = interpolate((rad_grid, chi_grid), cosas     ,method)
hits_interp    = interpolate((rad_grid, chi_grid), hits      ,method)

#wrapper for tan(phi) formalism
phi_interp_atan(rad,chi) = atan2(phi_interp_sin[rad,chi], phi_interp_cos[rad,chi])
