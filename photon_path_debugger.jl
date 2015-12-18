#Compare path of 2 photons


xpoint = 5.0
ypoint = -0.4

#xpoint = 0.01
#ypoint = -0.01
    
#ypoint = -1.0
#xpoint2 = 7.4374
xpoint2 = 5.0
#ypoint2 = -0.087
ypoint2 = -0.0017

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

println("timing bender...")
tic()
for i = 1:100
    rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                           beta, quad, wp, Rg)
    rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                           beta, quad, wp, Rg)
end
toc()

rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
println(pi/2-yns[end])
p1 = plot(rns, (pi/2-yns), "b-")#,xrange=[0.0, 0.0001])
rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p1 = oplot(rns, (pi/2-yns), "r--")
println(pi/2-yns[end])
    
rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p2 = plot(rns, zns, "b-")#, xrange=[0.0, 0.0001])
rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p2 = oplot(rns, zns, "r-")

    
rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p3 = plot(rns, ers, "b-", yrange=[0.0, 2.0e-5])
rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p3 = oplot(rns, ers, "r--")

rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                  beta, quad, wp, Rg)
p4 = plot(rns, lvs, "b-")
rns, yns, zns, ers, lvs, hit = bender3(xpoint2, ypoint2, sini,X, Osb,
                                  beta, quad, wp, Rg)
p4 = oplot(rns, lvs, "r-")

t = Table(4,1)
t[1,1] = p1
t[2,1] = p2
t[3,1] = p3
t[4,1] = p4    
    
