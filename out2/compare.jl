using Winston
using toolbox


#folder = "cadeau+morsink/"
folder = "f700/r12/"

#read file1
#fname1 = "f600_lamb_bb_R16.4_M1.4_rho1.csv"
fname1 = "f700_lamb_bb_R12.0_M1.4_rho10.csv"


da1 = readcsv(folder*fname1)
phase1 = da1[:,1]
flux1_kev_2 = da1[:,2] #Photon number flux at E (keV)
flux1_kev_6 = da1[:,3] #
flux1_kev_12 = da1[:,4] #
bNflux1 = da1[:,5] #Bolometric number flux ph/cm^2/s
bflux1 = da1[:,6] #Bolom flux keV/cm^2/s


#read file2
#fname2 = "resultcadfig3.dat"
#fname2 = "r16x10d90i45.dat"
#fname2 = "r12x10d90i45.dat"
fname2 = "r12x10d90i90.dat"


da2 = readdlm(folder*fname2)
phase2 = da2[:,1]
flux2_kev_2 = da2[:,2] #Photon number flux at E (keV)
flux2_kev_6 = da2[:,3] #
flux2_kev_12 = da2[:,4] #
bNflux2 = da2[:,5] #Bolometric number flux ph/cm^2/s
bflux2 = da2[:,6] #Bolom flux keV/cm^2/s



function comp_plot(phase1, flux1, phase2, flux2)
    
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
             ylabel="Flux (arb)"
             )
    p = oplot(phase2, flux2, "r.-")


    #interpolate abs. error
    Np = 128
    err = zeros(Np)
    for i = 1:Np
        iphase = phase2[i]
        val = toolbox.interp(phase1, flux1, iphase, method=:cubic)
        err[i] = (flux2[i] - val)/flux2[i]
    end

    pe = plot(phase2[1:Np], err,
              xrange = [0.0, 1.0],
              yrange = [-0.1, 0.1],
              xlabel = "Phase",
              ylabel = "Relative error"
              )

    pe = oplot([phase2[1], phase2[Np]], [0.0, 0.0], "k--")

    return p, pe
end

t = Table(4,3)

for i = 1:5
    if i == 1
        p1, pe1 = comp_plot(phase1, flux1_kev_2, phase2, flux2_kev_2)
        t[1,1] = p1
        t[2,1] = pe1
    elseif i == 2
        p2, pe2 = comp_plot(phase1, flux1_kev_6, phase2, flux2_kev_6)
        t[1,2] = p2
        t[2,2] = pe2
    elseif i == 3
        p3, pe3 = comp_plot(phase1, flux1_kev_12, phase2, flux2_kev_12)
        t[1,3] = p3
        t[2,3] = pe3
    elseif i == 4
        p4, pe4 = comp_plot(phase1, bNflux1, phase2, bNflux2)
        t[3,1] = p4
        t[4,1] = pe4
    elseif i == 5
        p5, pe5 = comp_plot(phase1, bflux1, phase2, bflux2)
        t[3,2] = p5
        t[4,2] = pe5
    #end
    elseif i == 6
        #interpolate ratio
        Np = 129
        ratio = zeros(Np)
        for i = 1:Np
            iphase = phase2[i]
            val = toolbox.interp(phase1, bflux1, iphase, method=:cubic)
            ratio[i] = (val)/bflux2[i]
        end
        pe6 = plot(phase2, ratio, "b-",
                   xlabel = "Phase",
                   ylabel = "ratio")
        pe6 = oplot([phase2[1], phase2[Np]], [1.0, 1.0], "k--")
        t[4,3] = pe6
    end
    
end

#t[1,1] = p1
#t[2,1] = p2

savefig(t, folder*fname1[1:end-3]*".eps")
