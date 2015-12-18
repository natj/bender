using Winston
using toolbox


#read file1
#fname1 = "f400_lamb_bb_R12.0_M1.6_rho30.csv"
#fname1 = "f400_lamb_bb_R12.0_M1.6_rho4.csv"
#fname1 = "f1_lamb_bb_R12.0_M1.6_rho1.csv"
#fname1 = "f1_lamb_bb_R12.0_M1.6_rho2.csv"
#fname1 = "f1_lamb_bb_R12.0_M1.6_rho4.csv"
fname1 = "f1_lamb_bb_R12.0_M1.6_rho30.csv"
#fname1 = "f1_lamb_bb_R12.0_M1.6_rho15.csv"

da1 = readcsv(fname1)

phase1 = da1[:,1]
bflux1 = da1[:,2]


#read file2
#fname2 = "nu400Hz_blackbody_rho30deg.dat"
#fname2 = "nu400Hz_blackbody_rho1deg.dat"
#fname2 = "nu1Hz_blackbody_rho1deg.dat"
fname2 = "nu1Hz_blackbody_rho30deg.dat"

da2 = readdlm(fname2)
phase2 = da2[:,1]
flux_2kev_2 = da2[:,2]
flux_2kev_6 = da2[:,3]
flux_2kev_12 = da2[:,4]
bNflux2 = da2[:,6]
bflux2 = da2[:,5]


#normalize
bflux1 = bflux1 ./ maximum(bflux1)
bflux2 = bflux2 ./ maximum(bflux2)


p1 = plot(phase1, bflux1, "b.-",
          xrange=[0.0, 1.0],
          yrange=[0.0, 1.1],
          xlabel="Phase",
          ylabel="Flux (arb)"
          )
p1 = oplot(phase2, bflux2, "r.-")


#interpolate abs. error
Np = 129
err = zeros(Np)
for i = 1:Np
    iphase = phase2[i]
    val = toolbox.interp(phase1, bflux1, iphase, method=:cubic)
    err[i] = (bflux2[i] - val)/bflux2[i]
end

p2 = plot(phase2[1:Np], err,
          xrange = [0.0, 1.0],
          yrange = [-0.1, 0.1],
          xlabel = "Phase",
          ylabel = "Relative error"
          )

p2 = oplot([phase2[1], phase2[Np]], [0.0, 0.0], "k--")

t = Table(2,1)
t[1,1] = p1
t[2,1] = p2

savefig(t, fname1[1:end-3]*".eps")
