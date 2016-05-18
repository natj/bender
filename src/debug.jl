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

println(pa)
println(pr)
println(pt)
println(pp)
