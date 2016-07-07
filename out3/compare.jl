using Winston
using toolbox

mass = 1.8 #mass
rad = 15 #radius
fs = 600 #frequency
incl = 30 #inclination

rho = 10 #spot radius
tx = 40 #spot colatitude
prof = "bb"  #beaming profile



folder = "my/"
#read file1
fname1 = "f$(round(Int,fs))p"*prof*"r$(round(Int,rad))m$(round(mass,1))d$(round(Int,tx))i$(round(Int,incl))x$(round(Int,rho)).csv"
println(fname1)

da1 = readcsv(folder*fname1)
phase1 = da1[:,1]
#flux1_kev_2 = da1[:,2] #Photon number flux at E (keV)
#flux1_kev_6 = da1[:,3] #
#flux1_kev_12 = da1[:,4] #
bNflux1 = da1[:,5] #Bolometric number flux ph/cm^2/s
bflux1 = da1[:,6] #Bolom flux keV/cm^2/s
flux1_kev_2 = da1[:,7] #Photon energy flux at E (keV)
flux1_kev_6 = da1[:,8] #
flux1_kev_12 = da1[:,9] #

#phase1 = phase1 .+ 0.03
#bflux1 = bflux1 .* 0.979
#flux1_kev_6 = flux1_kev_6 .* 7.5e10

#read file2
#f600r15m18x10d40i30_6keV_SD
fname2 = "f$(round(Int,fs))r$(round(Int,rad))m$(round(Int,10*mass))x$(round(Int,rho))d$(round(Int,tx))i$(round(Int,incl))_2keV_HT.csv"
da2 = readcsv(fname2)
phase2_kev_2 = da2[:,1]
flux2_kev_2 = da2[:,2] #Photon number flux at E (keV)

fname2 = "f$(round(Int,fs))r$(round(Int,rad))m$(round(Int,10*mass))x$(round(Int,rho))d$(round(Int,tx))i$(round(Int,incl))_6keV_HT.csv"
da2 = readcsv(fname2)
phase2_kev_6 = da2[:,1]
flux2_kev_6 = da2[:,2] #

fname2 = "f$(round(Int,fs))r$(round(Int,rad))m$(round(Int,10*mass))x$(round(Int,rho))d$(round(Int,tx))i$(round(Int,incl))_10keV_HT.csv"
da2 = readcsv(fname2)
phase2_kev_12 = da2[:,1]
flux2_kev_12 = da2[:,2] #


#wasted arrays
bNflux2 = da2[:,2] #Bolometric number flux ph/cm^2/s
bflux2 = da2[:,2] #Bolom flux keV/cm^2/s

#fmax1, imax1 = findmax(flux1_kev_2)
#fmax2, imax2 = findmax(flux2_kev_2)
#flux2_kev_2 = flux2_kev_2 .* fmax1/fmax2
#fmax1, imax1 = findmax(flux1_kev_12)
#fmax2, imax2 = findmax(flux2_kev_12)
#flux2_kev_12 = flux2_kev_12 .* fmax1/fmax2

#6kev
fmax1, imax1 = findmax(flux1_kev_6)
fmax2, imax2 = findmax(flux2_kev_6)
fmin1, imin1 = findmin(flux1_kev_6)
fmin2, imin2 = findmin(flux2_kev_6)
fshift = mean([fmax1/fmax2, fmin1/fmin2])
flux2_kev_6 = flux2_kev_6 .* fshift

#12kev
fmax1, imax1 = findmax(flux1_kev_12)
fmax2, imax2 = findmax(flux2_kev_12)
fmin1, imin1 = findmin(flux1_kev_12)
fmin2, imin2 = findmin(flux2_kev_12)
fshift = mean([fmax1/fmax2, fmin1/fmin2])
flux2_kev_12 = flux2_kev_12 .* fshift

#2kev
fmax1, imax1 = findmax(flux1_kev_2)
fmax2, imax2 = findmax(flux2_kev_2)
fmin1, imin1 = findmin(flux1_kev_2)
fmin2, imin2 = findmin(flux2_kev_2)
fshift = mean([fmax1/fmax2, fmin1/fmin2])
flux2_kev_2 = flux2_kev_2 .* fshift




#shift phase to match
#pshift = mean([ (phase1[imax1]-phase2_kev_6[imax2]),
#                          (phase1[imin1]-phase2_kev_6[imin2])])

pshift = -0.024
phase2_kev_2 = phase2_kev_2 .+ pshift
phase2_kev_6 = phase2_kev_6 .+ pshift
phase2_kev_12 = phase2_kev_12 .+ pshift


println("phase shift:",pshift)





function comp_plot(phase1, flux1, phase2, flux2, stitle)
    
    #normalize
    println()
    println("max  ref: $(maximum(flux2))")
    println("max comp: $(maximum(flux1))")
    println("   ratio: $(maximum(flux2)/maximum(flux1))")
    #flux1 = flux1 ./ maximum(flux1)
    #flux2 = flux2 ./ maximum(flux2)

    p = plot(phase1, flux1, "b.-",
             xrange=[0.0, 1.0],
             #yrange=[0.0, 1.0],
             xlabel="Phase",
             #ylabel="Flux (arb)"
             ylabel=stitle
             )
    p = oplot(phase2, flux2, "r.-")


    #interpolate abs. error
    Np = length(phase2)
    err = zeros(Np)
    for i = 1:Np
        iphase = phase2[i]
        val = toolbox.interp(phase1, flux1, iphase, method=:cubic)
        err[i] = (flux2[i] - val)/flux2[i]
    end

    pe = plot(phase2[1:Np], 100.* err,
              xrange = [0.0, 1.0],
              #yrange = [-0.2, 0.2],
              #yrange = [-0.1, 0.1],
              xlabel = "Phase",
              ylabel = "Relative error (%)"
              )

    pe = oplot([phase2[1], phase2[Np]], [0.0, 0.0], "k--")

    return p, pe
end

t = Table(4,3)


# for i = 1:5
#     if i == 2
#         p2, pe2 = comp_plot(phase1, flux1_kev_6, phase2, flux2_kev_6, "N (6 keV) [ph/cm^2/s/keV]")
#         t[1,2] = p2
#         t[2,2] = pe2
#     end
# end
for i = 1:5
    if i == 1
        p1, pe1 = comp_plot(phase1, flux1_kev_2, phase2_kev_2, flux2_kev_2, "N (2 keV) [ph/cm^2/s/keV]")
        t[1,1] = p1
        t[2,1] = pe1
    elseif i == 2
        p2, pe2 = comp_plot(phase1, flux1_kev_6, phase2_kev_6, flux2_kev_6, "N (6 keV) [ph/cm^2/s/keV]")
        t[1,2] = p2
        t[2,2] = pe2
    elseif i == 3
        p3, pe3 = comp_plot(phase1, flux1_kev_12, phase2_kev_12, flux2_kev_12, "N (12 keV) [ph/cm^2/s/keV]")
        t[1,3] = p3
        t[2,3] = pe3
    #elseif i == 4
        #p4, pe4 = comp_plot(phase1, bNflux1, phase2, bNflux2, "Bolometric N [ph/cm^2/s]")
        #t[3,1] = p4
        #t[4,1] = pe4
    #elseif i == 5
        #p5, pe5 = comp_plot(phase1, bflux1, phase2, bflux2, "Bolometric E [erg/cm^2/s]")
        #t[3,2] = p5
        #t[4,2] = pe5
    end
end

#t[1,1] = p1
#t[2,1] = p2

savefig(t, folder*fname1[1:end-3]*".eps")
