#Compare path of 2 photons
include("bender.jl")

xpoint = 0.1
ypoint = 7.4

#xpoint = 0.01
#ypoint = -0.01
    
#ypoint = -1.0
#xpoint2 = 7.4374
xpoint2 = 0.1
#ypoint2 = -0.087
ypoint2 = -6.8

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

# println("timing bender...")
# tic()
# for i = 1:100
#     rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
#                                            beta, quad, wp, Rg)
#     rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
#                                            beta, quad, wp, Rg)
# end
# toc()

tcor = G*M/c^3

rns, yns, zns, tns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                       beta, quad, wp, Rg)
println("x y end $(rns[end]) $(tns[end]*tcor)")
println(pi/2-yns[end])

p1 = plot(rns, (pi/2-yns), "b-")#,xrange=[0.2, 0.21])
p2 = plot(rns, zns, "b-")#, xrange=[0.0, 0.0001])
p3 = plot(rns, ers, "b-", yrange=[0.0, 2.0e-5])
p4 = plot(rns, lvs, "b-")


##########
rns, yns, zns, tns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                       beta, quad, wp, Rg)
println("x y end $(rns[end]) $(tns[end]*tcor)")
println(pi/2-yns[end])

p1 = plot(p1, rns, (pi/2-yns), "r--")
p2 = plot(p2, rns, zns, "r-")
p3 = plot(p3, rns, ers, "r--")
p4 = plot(p4, rns, lvs, "r-")




t = Table(4,1)
t[1,1] = p1
t[2,1] = p2
t[3,1] = p3
t[4,1] = p4    
    
