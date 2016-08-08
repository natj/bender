using Cubature

#Compute spot on the NS
include("plot2d.jl")

#Load photon bender
#include("bender.jl")

#Compute raw image
#include("rtrace.jl")

#Interpolate from raw image and compute radiation processes
#include("radiation.jl")

rho = deg2rad(30.0)
colat = deg2rad(50.0)

interp = true
exact_edges = false


########
function spot(t, phi, theta;
              stheta = deg2rad(50.0), #spot latitude
              delta = deg2rad(30.0) #spot half-opening angle
              )

    #Circular spot
    d = great_circle_dist(0.0, phi, stheta, theta)
    if abs(d) < delta
        return true
    end

    #secondary spot
    #d = great_circle_dist(0.0, phi+pi, pi-stheta, theta)
    #if abs(d) < delta
    #    return true
    #end

    
    return false
end


img4 = zeros(Ny_dense, Nx_dense) #debug array

#Spot image frame size
N_frame = 20

N_frame_chi = 500
N_frame_rad = 100

#Beaming function for the radiation
#Ir(cosa) = 1.0 #isotropic beaming
#Ir(cosa) = cosa

#Time parameters
Nt = 64

times = collect(linspace(0, 1/fs, Nt))
tbin = abs(times[2] - times[1])/2.0 
phase = collect(times .* fs)

spot_flux = zeros(Nt)
spot_flux2 = zeros(Nt)

sdelta = zeros(Nt)
sdelta2 = zeros(Nt)
sdelta3 = zeros(Nt)

sfluxE = zeros(Nt, 3)
sfluxNE = zeros(Nt, 3)

sfluxB = zeros(Nt)
sfluxNB = zeros(Nt)

iters = zeros(1)

function dF(xx, yy)

    dfluxNB = 0.0
    dfluxB = 0.0
    dfluxE = zeros(3)
    dfluxNE = zeros(3)


    x = xx[1]
    y = xx[2]
    #println("dF: x:$x y:$y $imgscale $t")

    rad = hypot(x,y)
    chi = mod2pi(pi/2 - atan2(y,x))

    if interp
        phi = phi_interp_atan(rad,chi)
        theta = theta_interp[rad,chi]
        Xob = Xs_interp[rad,chi]
        time = time_interp[rad,chi]
        cosa = cosa_interp[rad,chi]
        dtau = dtau_interp[rad,chi]
        #hit = hits_interp[rad,chi] #test if we hit the surface
    
        #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
    else
        time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                    X, Osb, beta, quad, wp, Rg)
        time -= time0
        #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
        #println()
    end

    if rad <= edge_interp(chi)
        dt = time*G*M/c^3
        phi = phi - (t - dt)*fs*2*pi
        phi = mod2pi(phi)
        inside = spot(0.0, phi, theta,
                      stheta = colat,
                      delta = rho
                      )
        if inside
    
            if interp
                delta = delta_interp[rad,chi]
                EEd = reds_interp[rad,chi]
            #println(EEd, " ", delta)
            else
    
                # update subimage corners; they might no be up-to-date
                # because we trace exactly here instead of interpolating.
                #old_subframe[1] = frame_y2 < y ? y : frame_y2 #top #max
                #old_subframe[2] = frame_y1 > y ? y : frame_y1 #bottom #min
                #old_subframe[3] = frame_x1 > x ? x : frame_x1 #left min
                #old_subframe[4] = frame_x2 < x ? x : frame_x2 #right max
    
                
                EEd, delta, dtau = radiation(rad, chi,
                                             phi, theta, cosa,
                                             X, Xob, Osb, sini)
            end
            #println(EEd, " ", delta)
            #println()
    
            if exact_edges && delta < 0.98
                #println("exact edge for $rad")
                time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                            X, Osb, beta, quad, wp, Rg)
                time -= time0
                dt = time*G*M/c^3
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho)
                if inside && hit
                    EEd, delta, dtau = radiation(rad, chi,
                                                 phi, theta, cosa,
                                                 X, Xob, Osb, sini)
    
                else
                    EEd = 1.0
                    delta = 1.0
                    dtau = 0.0
                end
            end
    
            dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)

    
            dfluxB *= dtau
            dfluxNB *= dtau
            for ie = 1:3
               dfluxE[ie] *= dtau 
               dfluxNE[ie] *= dtau
            end
    
        end #inside spot
    end#hiti


    yy[:] = [dfluxE[1],dfluxE[2],dfluxE[3],
            dfluxNE[1],dfluxNE[2],dfluxNE[3],
            dfluxNB,
            dfluxB
            ]

    iters[1] += 1

    return 1
end




#Cartesian integration
#########
old_subframe = [y_grid_d[1],
                y_grid_d[end],
                x_grid_d[end],
                x_grid_d[1]
                ]

t = 0.0
frame_dxdy = 1.0


for k = 1:Nt
#for k = 1:1
    img4[:,:] = 0.0
    
    t = times[k]
    println()
    println("t: $t k: $k")
    
    #set whole image as starting frame
    frame_y2 = y_grid_d[1]
    frame_y1 = y_grid_d[end]
    frame_x1 = x_grid_d[end]
    frame_x2 = x_grid_d[1]

    located_spot = false
    
    for j = y1s:y2s
        y = y_grid_d[j]

        #frame_left = false
        
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))

            if rad <= edge_interp(chi)

                #trace back to star
                phi = phi_interp_atan(rad,chi)
                theta = theta_interp[rad,chi]
                time = time_interp[rad,chi]
                dtau = dtau_interp[rad,chi]
                
                #rotate star
                dt = time*G*M/c^3 #time shift
                #dt = 0.0
                
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
      
                #img4[j,i] = -4.0*painter(phi, theta)/2.0
                img4[j,i] = painter(phi, theta)/2.0
                
                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho
                              )
             
                if inside 
                    located_spot = true
                    #track down spot edges

                    #println("x = $x y = $y")
                    frame_y2 = frame_y2 < y ? y : frame_y2 #top #max
                    frame_y1 = frame_y1 > y ? y : frame_y1 #bottom #min
                    frame_x1 = frame_x1 > x ? x : frame_x1 #left min
                    frame_x2 = frame_x2 < x ? x : frame_x2 #right max

                    #println(frame_x1)
                    #println(frame_x2)
                    #println(frame_y1)
                    #println(frame_y2)
                    #println()
                     
                    #Time shifts for differnt parts
                    time = time_interp[rad,chi]
                    cosa = cosa_interp[rad,chi]
    
                    delta = delta_interp[rad,chi]
                    EEd = reds_interp[rad,chi]
                        
                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)
                        
                    img4[j,i] += dfluxB * dxdy * imgscale*1.0e7
                    
                end #if inside            

            
            end #hiti
            #end #XXX debug
        end# for x
    end#for y

    if !located_spot
        println("hidden spot")
        frame_y2 = old_subframe[1]
        frame_y1 = old_subframe[2]
        frame_x1 = old_subframe[3]
        frame_x2 = old_subframe[4]
    end

    #expand image a bit
    #########
    frame_expansion_x = abs(x_grid_d[8] - x_grid_d[1])
    frame_expansion_y = abs(y_grid_d[8] - y_grid_d[1])
    frame_y2 += frame_expansion_y
    frame_y1 -= frame_expansion_y
    frame_x1 -= frame_expansion_x
    frame_x2 += frame_expansion_x

    #Exapand image a bit keeping the aspect ratio
    ##########
    frame_expansion_x2 = abs(frame_x2 - frame_x1)
    frame_expansion_y2 = abs(frame_y2 - frame_y1)

    frame_xs = abs(frame_x2 - frame_x1)/N_frame
    frame_ys = abs(frame_y2 - frame_y1)/N_frame


    println("x1: $frame_x1  x2: $frame_x2  y1: $frame_y1 y2: $frame_y2")  
    println("x = $frame_xs y = $frame_ys")

    #pick smaller
    if frame_xs > frame_ys
        Nx_frame = N_frame
        Ny_frame = max(round(Int, (frame_ys*N_frame/frame_xs)), 2)
    else
        Ny_frame = N_frame
        Nx_frame = max(round(Int, (frame_xs*N_frame/frame_ys)), 2)
    end

    println("Nx= $Nx_frame Ny = $Ny_frame")


    frame_xgrid = collect(linspace(frame_x1, frame_x2, Nx_frame))
    frame_ygrid = collect(linspace(frame_y1, frame_y2, Ny_frame))
    frame_dxx = 1.0*(frame_xgrid[2] - frame_xgrid[1])
    frame_dyy = 1.0*(frame_ygrid[2] - frame_ygrid[1])
    frame_dxdy = frame_dxx*frame_dyy #*X^2

    #Locate spot edges on the old grid
    ##########
    
    #Plot large image with bounding box for the spot
    p10a = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 2, "Blues")
    Winston.add(p10a, Curve([frame_x1, frame_x2, frame_x2, frame_x1, frame_x1],
                           [frame_y1, frame_y1, frame_y2, frame_y2, frame_y1],
                           linestyle="solid"))
    #add time stamp
    xs = x_grid_d[1] + 0.84*(x_grid_d[end] - x_grid_d[1])
    ys = y_grid_d[1] + 0.93*(y_grid_d[end] - y_grid_d[1])
    Winston.add(p10a, Winston.DataLabel(xs, ys, "$(k) ($(round(times[k]*fs,3)))"))
    #display(p10a)

    #println("dx = $(frame_dxx) dy = $(frame_dyy)")



    img5 = zeros(Ny_frame, Nx_frame) #debug array

    Ndelta = 0.0

    old_subframe = [frame_y2,
                    frame_y1,
                    frame_x1,
                    frame_x2
                    ]

    #Integration
    tic()
    iters[1] = 0
    for j = 1:Ny_frame
        y = frame_ygrid[j]
        for i = 1:Nx_frame
            x = frame_xgrid[i]

            dFs = zeros(8)
            istat = dF([x,y], dFs) 
            dFs *= imgscale * frame_dxdy

            sfluxE[k, 1] += dFs[1]
            sfluxE[k, 2] += dFs[2]
            sfluxE[k, 3] += dFs[3]

            sfluxNE[k, 1] += dFs[4]
            sfluxNE[k, 2] += dFs[5]
            sfluxNE[k, 3] += dFs[6]

            sfluxNB[k] += dFs[7]
            sfluxB[k] += dFs[8]
        end #x
    end#y
    println("Riemann integral iterations: $(iters[1])")
    toc()

    #println("bol flux: $(sfluxB[k]) | num flux $(sfluxNB[k])")
    println("num flux $(sfluxNB[k])")

    tic()
    iters[1] = 0
    (vals, errs) = pcubature(8,
                           dF,
                           [frame_x1, frame_y1],
                           [frame_x2, frame_y2],
                           reltol = 1.0e-3,
                           abstol = 0.0,
                           maxevals=0)
    println("Cubature iterations: $(iters[1])")
    toc()

    vals *= imgscale
    errs = errs ./ vals #transform into relative err
    errs *= imgscale

    val = vals[7]
    err = errs[7]
    merr = (sfluxNB[k] - val)/val

    sfluxE[k, 1] = vals[1]
    sfluxE[k, 2] = vals[2]
    sfluxE[k, 3] = vals[3]

    sfluxNE[k, 1] = vals[4]
    sfluxNE[k, 2] = vals[5]
    sfluxNE[k, 3] = vals[6]

    sfluxNB[k] = vals[7]
    sfluxB[k] = vals[8]

    println("num flux $val $err $merr")
    p10b = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 2, "Blues")

    #add time stamp
    xs = frame_xgrid[1] + 0.84*(frame_xgrid[end]-frame_xgrid[1])
    ys = frame_ygrid[1] + 0.93*(frame_ygrid[end]-frame_ygrid[1])
    Winston.add(p10b, Winston.DataLabel(xs, ys, "$(k)"))
    #display(p10)

    p10c = plot(phase, sfluxNB[:], "k-")
    p10c = oplot([phase[k]], [sfluxNB[k]], "ko")

    #doppler factor
    sdelta[k] = sdelta3[k]/Ndelta
    sdelta2[k] = sdelta3[k]/Ndelta
    sdelta3[k] = sdelta3[k]/Ndelta

    sdelta[k] = sfluxNB[k]/val

    p10d = plot(phase, sdelta, "b-",
                xrange=[0.0, 1.0],
                yrange=[0.8, 1.2])




    #Plot large image with more fine-tuned aesthetics
    #mkdir("movie/")
    #p10aa = plot2d(img4, x_grid_d, y_grid_d, 0, -8, 8, "RdBu")
    #p10aa = plot2d(clamp(img4,0.0,1.9), x_grid_d, y_grid_d, 0, 0, 2, "Blues")
    #Winston.setattr(p10aa.frame, 
    #                draw_spine=false, 
    #                draw_ticks=false,
    #                draw_ticklabels=false)
    #savefig(p10aa, "movie/bender_$(lpad(k-1,3,"0")).png")

    tt1 = Table(1,2)
    tt1[1,1] = p10a
    tt1[1,2] = p10b

    tt = Table(3,1)
    tt[1,1] = tt1
    tt[2,1] = p10c
    tt[3,1] = p10d
    display(tt)

    #readline(STDIN)
end#for t



#write to file
opath = "out_cub/"
#opath = "out2/"
#opath = "out2/cadeau+morsink/"
#opath = "out2/f$(round(Int,fs))/r$(round(Int,R/1e5))n/"
#opath = "out3/HT/"
#opath = "out4/"

mkpath(opath)

fname = "f$(round(Int,fs))pbbr$(round(Int,R/1e5))m$(round(M/Msun,1))d$(round(Int,rad2deg(colat)))i$(int((rad2deg(incl))))x$(round(Int,rad2deg(rho))).csv"
#fname = "f$(round(Int,fs))phopfr$(round(Int,R/1e5))m$(round(M/Msun,1))d$(round(Int,rad2deg(colat)))i$(int((rad2deg(incl))))x$(round(Int,rad2deg(rho))).csv"


wmatr = zeros(Nt, 9)
wmatr[:,1] = phase
wmatr[:,2] = sfluxNE[:, 1] #2 kev
wmatr[:,3] = sfluxNE[:, 2] #6 kev
wmatr[:,4] = sfluxNE[:, 3] #12 kev
wmatr[:,5] = sfluxNB #bol number flux
wmatr[:,6] = sfluxB #bol energy flux
wmatr[:,7] = sfluxE[:, 1] #2 kev
wmatr[:,8] = sfluxE[:, 2] #6 kev
wmatr[:,9] = sfluxE[:, 3] #12 kev

writecsv(opath*fname, wmatr)



