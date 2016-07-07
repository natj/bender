include("../bender.jl")

rrs = 0.2
yis = 0.4

nu2   = beta/3.0 - quad*0.5*(3*cos(yis)^2-1)
B2    = beta
zeta2 = beta*(3*0.5*(3*cos(yis)^2-1)/4-1/3)

x = 2.5
y = 3.0

pa        = ptim(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
pr, rturn = prad(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
pt, tturn = pthe(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
pp        = pphi(x, y, sini, rrs, nu2, B2, zeta2, wp, yis, Rg)
dr        = -Rg/rrs^2 #dx to dr

println("tim ", pa)
println("rad ", pr)
println("the ", pt)
println("phi ", pp)


t0, phi, theta, Xob, hit, cosa = bender3(x, y, sini, X, Osb, beta, quad, wp, Rg)
#rns, yns, zns, tns, ers, lvs, hit  = bender3(x, y, sini, X, Osb, beta, quad, wp, Rg)


println("timec ", t0)
println("rad   ", Xob)
println("the   ", theta)
println("phi   ", phi)
println("cosa  ", cosa)


