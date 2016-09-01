using Winston
using Cubature
using Interpolations

#Rayleigh (Thomson) scattering
psi(mu) = (3/16)*(3-mu^2)

#for polarization
psi1(mu) = (3/4)*(1-mu^2)
psi2(mu) = (3/8)*(1-mu^2)



################################
method = Gridded(Linear())

#first guess
HnA(mu) = 1.0 .+ 2.3.*mu .- 0.3.*mu.^2
Hnm1(mu) = HnA(mu)

mus = collect(linspace(0.01, 1.0, 1000))

Hns = zeros(length(mus))
Hns[:] = HnA(mus)

rvalt = zeros(length(mus))
amoms = zeros(3)

p = plot(xrange=[0,1],
         yrange=[0.35, 1.3],
         xlabel="mu",
         ylabel="H(mu)/2a")

for iter = 1:10
    
    #println()
    println("iter $iter")

    #make array of H values into function by interpolation
    #println(Hnm1(mus))
    Hn_interp = interpolate((mus,), Hns, method)
    Hnm1(eta) = Hn_interp[eta]

    #calculate moments of H
    for mom in [0,1,2]
        (val0, err0) = pcubature(x -> Hnm1(x[1])*x[1]^mom,
                                0,1,reltol=1.0e-8,maxevals=0)
        amoms[mom+1] = val0
    end
    println("Moments: a0: $(amoms[1]) a1: $(amoms[2]) a2: $(amoms[3])")
    

    for i = 1:length(mus)
        
        #loop over mu values
        mu = mus[i]

        #first term
        integ1(x) = psi(x[1])
        (val1, err) = pcubature(integ1,0,1,reltol=1.0e-8,maxevals=0)

        #second term
        (val2, err2) = pcubature(x -> begin x[1]*psi(x[1])/(mu+x[1])*Hnm1(x[1]); end,
                                0,1,reltol=1.0e-8,maxevals=0)

        rvalt[i] = sqrt(1.0 - 2.0*val1) + val2
        println("T1: $val1 T2: $val2 I: $(rvalt[i])")
    end

    #Hns[:] = 1.0 ./ rvalt
    Hns[:] = 0.5*((1.0 ./ rvalt) .+ Hns)

    p = oplot(mus, Hns./2./amoms[2], "k-")
    display(p)
    #readline(STDIN)
end    
    

#calculate moments of H
function Hn_moments(Hnf)
    moms = zeros(3)
    for mom in [0,1,2]
        (val0, err0) = pcubature(x -> Hnf(x[1])*x[1]^mom,
                                0,1,reltol=1.0e-8,maxevals=0)
        moms[mom+1] = val0
    end
    return moms
end

#Solution of the integral equation
Hne(mu) = Hnm1(mu)
momse = Hn_moments(Hne)
refval = Hne(mus)./(2*momse[2])
p = plot(mus, Hne(mus)./(2*momse[2]), "k-")
erre = (Hne(mus)./(2*momse[2])) ./ refval -1

#1st order approx
Hn1(mu) = 1 + 2.06.*mu
moms1 = Hn_moments(Hn1)
p = oplot(mus, Hn1(mus)./(2*moms1[2]), "m--")
err1 = (Hn1(mus)./(2*moms1[2])) ./ refval -1

#2nd order approx
Hn2(mu) = 1 .+ 2.3.*mu .- 0.3*mu.^2
moms2 = Hn_moments(Hn2)
p = oplot(mus, Hn2(mus)./(2*moms2[2]), "g--")
err2 = (Hn2(mus)./(2*moms2[2])) ./ refval -1

#FL approx 
HnAFL(mu) = 0.42822 + 0.92236.*mu - 0.085751.*mu.^2
p = oplot(mus, HnAFL(mus), "b--")
errFL = (HnAFL(mus)) ./ refval -1

#p2 = plot(xrange=[0,1],
#          yrange=[-0.05, 0.05],
#          xlabel="mu",
#          ylabel="ratio-1")
#p2 = oplot(mus, erre, "k-")
#p2 = oplot(mus, err1, "m-")
#p2 = oplot(mus, err2, "g-")
#p2 = oplot(mus, errFL, "b-")
#display(p2)



#Final beaming function
#Beaming function
Beam(mu) = 1.0 #Lambertian
#Beam(mu) = 0.42822 + 0.92236*mu - 0.085751*mu^2 #approx Hopf
#Beam(mu) = (1.0 + 2.3*mu - 0.3*mu^2)/(2*1.194)
#Beam(mu) = Hne(mu)./(2*momse[2])

