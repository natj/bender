#make synthetic data out of bender output
using Winston

#fpath = "out3/HT/f400pbbr13m1.5d50i60x10.csv"
fpath = "out2/f600/f600pbbr15m1.6d50i60x30_obl_gamma_test2.csv"

rmatr = readcsv(fpath)
phase    = rmatr[:,1]
fluxNE2  = rmatr[:,2]
fluxNE6  = rmatr[:,3]
fluxNE12 = rmatr[:,4]

fluxNB   = rmatr[:,5]
fluxB    = rmatr[:,6]
fluxE2   = rmatr[:,7]
fluxE6   = rmatr[:,8]
fluxE12  = rmatr[:,9]


#assume channels
#1.0 - 3.0   #dE = 1.0 keV # Emean = 2.0
#3.0 - 9.0   #dE = 3.0 keV # Emean = 6.0
#9.0 - 15.0  #dE = 3.0 keV # Emean = 12.0
Emeans = [2.0, 6.0, 12.0]
dEs = [2.0, 6.0, 6.0]

#compute observed seconds
fs = 400
dur = 1.0/fs

#0 -> 0.1 -> 0.2 -> 0.3
# 0.05   0.05

#number of photons is N_E x dE x dt
# and so we observe ph/cm^2

obsF = zeros(length(phase), 3)
expo = zeros(length(phase))

tb = 0.0
pb = 0.0

dphase = (phase[2] - phase[1])/2
for t in 1:length(phase)
    
    d1 = dphase
    
    #if !(t == 1 || t == length(phase))
        d1 *= 2.0
    #end

    pb += d1
    tb += d1*dur
    
    println(tb, " ", pb)

    obsF[t, 1] = d1*dur*dEs[1]*fluxNE2[t]
    obsF[t, 2] = d1*dur*dEs[2]*fluxNE6[t]
    obsF[t, 3] = d1*dur*dEs[3]*fluxNE12[t]
    expo[t] = d1*dur
end

dErr = 1.0e-7
yerr1 = dErr .+ abs(0.75*dErr*randn(length(phase)))
yerr2 = dErr .+ abs(0.75*dErr*randn(length(phase)))
yerr3 = dErr .+ abs(0.75*dErr*randn(length(phase)))



obsF[:,1] += 0.1*dErr*randn(length(phase))
obsF[:,2] += 0.1*dErr*randn(length(phase))
obsF[:,3] += 0.1*dErr*randn(length(phase))


# yerr1 = yerr1/maximum(obsF[:,1])
# obsF[:,1] = obsF[:,1]/maximum(obsF[:,1])
# yerr2 = yerr2/maximum(obsF[:,2])
# obsF[:,2] = obsF[:,2]/maximum(obsF[:,2])
# yerr3 = yerr3/maximum(obsF[:,3])
# obsF[:,3] = obsF[:,3]/maximum(obsF[:,3])



p = plot(phase, obsF[:,1], "b-")
p = errorbar(p, phase, obsF[:,1], yerr=yerr1, color="blue")
oplot(phase, obsF[:,2], "r-")
p = errorbar(p, phase, obsF[:,2], yerr=yerr2, color="red")
oplot(phase, obsF[:,3], "g-")
p = errorbar(p, phase, obsF[:,3], yerr=yerr3, color="green")


f = open("synt_data.dat", "w")

print(f,"# Data is ph/cm^2
#
# one spot
# D = 10 kpc
# nu = 600 Hz
# m 1.0 - 2.0 Msun
# r 8.0 - 18.0 km
# rho = 5 - 40 deg
# incl = 0 - 90 deg
# theta_s = 0 - 90 deg
#
#
# columns
# 1 Phase
# 2 exposure (s) 
# 3 number of photons in channel 1.0 - 3.0 keV  (dE = 2.0 keV; Emid = 2.0 keV)
# 4                              3.0 - 9.0 keV  (dE = 6.0 keV; Emid = 6.0 keV)
# 5                              9.0 - 15.0 keV (dE = 6.0 keV; Emid = 12.0 keV)
# 6 error in channel             1.0 - 3.0 keV
# 7                              3.0 - 9.0 keV
# 8                              9.0 - 15.0 keV
#
")


for i = 1:length(phase)
    @printf(f,"%5.3e %9.6e %9.6e %9.6e %9.6e %9.6e %9.6e %9.6e\n",
            phase[i],
            expo[i],
            obsF[i, 1],
            obsF[i, 2],
            obsF[i, 3],
            yerr1[i],
            yerr2[i],
            yerr3[i]
            )
end

close(f)
