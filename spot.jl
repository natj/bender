#Compute spot on the NS

include("plot2d.jl")

#Load photon bender
#include("bender.jl")

#Compute raw image
#include("rtrace.jl")

#Interpolate from raw image and compute radiation processes
#include("radiation.jl")

rho = deg2rad(1.0)
colat = deg2rad(50.0)


interp = true

########
function spot(t, phi, theta;
              stheta = deg2rad(50.0), #spot latitude
              delta = deg2rad(30.0) #spot half-opening angle
              )

    #Vincenty's formula
    d = great_circle_dist(0.0, phi, stheta, theta)
    if abs(d) < delta
        return true
    end

    #secondary spot
    # d = great_circle_dist(0.0, phi+pi, pi-stheta, theta)
    # if abs(d) < delta
    #     return true
    # end

    
    return false
end


#photon time lag function
# returns new time-array index taking into account time delay effect
function time_lag(t, k, times, Nt, tbin, phi, theta)

    #exact from raytracing
    dt = t*G*M/c^3
    
    #approximative
    #cosi = sqrt(1-sini^2)
    #cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
    #y = 1 - cospsi

    #dt = y*R/c
    #dt = y*(1.0 + (U*y/8.0)*(1+ y*(1.0/3.0 - U/14.0)))*R/c
    #println("dt: $dt  dt0: $dt0  dt2: $dt2 ")
    #println(" $(dt01/dt) ")

    #get new timebin index
    kd = 1
    while dt >= (2kd - 1)*tbin
        kd += 1
    end
    kd -= 1

    kindx = k + kd
    if kindx > Nt; kindx -= Nt; end
    #while kindx > Nt
    #    kindx -= Nt
    #end
    
    #println("k: $k kd: $kd kindx: $kindx")
    #println()
    
    return kindx
end


img4 = zeros(Ny_dense, Nx_dense) #debug array

#Spot image frame size
N_frame = 500


#Beaming function for the radiation
#Ir(cosa) = 1.0 #isotropic beaming
#Ir(cosa) = cosa

#Time parameters
Nt = 128

times = collect(linspace(0, 1/fs, Nt))
tbin = abs(times[2] - times[1])/2.0 
phase = collect(times .* fs)

spot_flux = zeros(Nt)
spot_flux2 = zeros(Nt)

sdelta = zeros(Nt)
sdelta2 = zeros(Nt)

sfluxE = zeros(Nt, 3)
sfluxNE = zeros(Nt, 3)

sfluxB = zeros(Nt)
sfluxNB = zeros(Nt)

#Polar integration

#create thick radius grid
Nrad_frame = 1000
rad_diffs = 1 ./ exp(linspace(0.0, 1.5, Nrad_frame-1).^2)
rad_grid_d = rmax * cumsum(rad_diffs) / sum(rad_diffs)
unshift!(rad_grid_d, 0.0)

tic()
#for k = 1:Nt
for k = 2:1
    t = times[k]

    radmin = 0.0
    radmax = 10.0
    chimin = 0.0
    chimax = 2pi
    
    for i = 2:Nchi-1
        chi = chi_grid[i]

        #avoid double integrating angle edges
        if chi < 0; continue; end
        if chi > 2pi; continue; end
        
        #if mod(i,10) == 0; println("chi: $(round(chi/2pi,2))"); end
    
        for j = 2:Nrad-1
            rad = rad_grid[j]

            if hits[j, i] < 1; break; end

    
            #Ray traced photons
            ####
            phi = Phis[j, i]
            theta = Thetas[j, i]

            #rotate star
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)

            inside = spot(0.0, phi, theta,
                          stheta = colat,
                          delta = rho
                          )

            
            if inside #|| inside2
                #println("inside")
                radmin = rad > radmin ? rad : radmin
                radmax = rad < radmax ? rad : radmax
                                 
                #time = Times[j, i]
                #Xob = Xs[j, i]
                #cosa = cosas[j, i]
                #dF = dFlux[j, i]
                #dE = Reds[j, i]

                #kd = time_lag(time, k, times, Nt, tbin, phi, theta)

                #drdchi = abs(rad_grid[j+1] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
                #spot_flux2[kd] += dF * drdchi
            end#inside
        end#rad
    end#chi


    #Integrate in thick interpolated grid
    println("radmin: $radmin radmax: $radmax")

       
    for i = 2:Nchi-1
        chi = chi_grid[i]

        #avoid double integrating angle edges
        if chi < 0; continue; end
        if chi > 2pi; continue; end
        
        for j = 2:Nrad_frame-1
        #j = 2
        #while j <= Nrad_frame-1    
            rad = rad_grid_d[j]

            #if rad < radmin; continue; end
            #if rad > radmax; j = Nrad_frame; end
            
            
            if hits_interp[rad, chi] < 1; break; end

            phi = phi_interp_atan(rad, chi)
            theta = theta_interp[rad, chi]
            
            #rotate star
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)
            
            inside = spot(0.0, phi, theta,
                          stheta = colat,
                          delta = rho
                          )

            #println("spot")
            
            if inside
                #println("inside")
                time = time_interp[rad, chi]
                Xob = Xs_interp[rad, chi]
                cosa = cosa_interp[rad, chi]
                dF = flux_interp[rad, chi]
                dE = reds_interp[rad, chi]

                kd = time_lag(time, k, times, Nt, tbin, phi, theta)

                drdchi = abs(rad_grid_d[j+1] - rad_grid_d[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid_d[j]
                spot_flux2[kd] += dF * drdchi
            end#inside

            #j += 1
        end
    end

    
    println("time = $t")
    p10polar = plot(phase, spot_flux2, "k-")
    p10polar = oplot([phase[k]], [spot_flux2[k]], "ko")
    display(p10polar)
    
end#time
toc()


#Cartesian integration
#########
old_subframe = [y_grid_d[1],
                y_grid_d[end],
                x_grid_d[end],
                x_grid_d[1]
                ]


tic()

for k = 1:Nt
#for k = 55:Nt
#for k = 12:23
#for k = 40:40
#for k = 27:26
#for k = 80:80
#for k = 24:38
#for k = 20:45
        
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
    #inf_small = true
    
    for j = y1s:y2s
        y = y_grid_d[j]

        #frame_left = false
        
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))

            #hit = hits_interp[rad,chi] #test if we hit the surface
            #hiti = round(Int, hit)
            #hiti = rad < rmax ? hiti : 0
            
            #if (rad > 5) && (pi/3.5 < chi < pi/2)
            #if hiti > 0
            if rad <= edge_interp(chi)
            
            #trace back to star
            phi = phi_interp_atan(rad,chi)
            theta = theta_interp[rad,chi]
            time = time_interp[rad,chi]
            
            
            #rotate star
            dt = time*G*M/c^3 #time shift
            #dt = 0.0
            
            phi = phi - (t -dt)*fs*2*pi
            phi = mod2pi(phi)
      
            #img4[j,i] = -4.0*painter(phi, theta)/2.0
            img4[j,i] = painter(phi, theta)/2.0
            
            inside = spot(0.0, phi, theta,
                          stheta = colat,
                          delta = rho
                          )
            
            if inside
                located_spot = true
            #if inside && inf_small
            #    inf_small = false
                
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
                kd = time_lag(time, k, times, Nt, tbin, phi, theta)
    
                #Xob = Xs_interp[y,x] 
                #cosa = cosa_interp[y,x]
                #dF, dE = radiation(Ir,
                #                   x,y,
                #                   phi, theta, cosa,
                #                   X, Xob, Osb, sini, img3[j,i])
                #dF = flux_interp[y,x]
                #dE = reds_interp[y,x]
                delta = delta_interp[rad,chi]
                EEd = reds_interp[rad,chi]
                    
                dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)
                    
                #img4[j,i] = painter(phi, theta)
                    
                #img5[j,i] += 1.0e9*dF * frame_dxdy
                #println(dF)

                #println(dfluxB)
                #img4[j,i] = abs(img4[j,i]) + dfluxB * dxdy * imgscale*1.0e7
                img4[j,i] += dfluxB * dxdy * imgscale*1.0e7
                #img4[j,i] += EEd * dxdy * imgscale /1.0e5
                
                #println(img4[j,i])
                #dF = flux_interp[rad,chi]
                #dE = reds_interp[rad,chi]
                                                
                #img4[j,i] = painter(phi, theta)
                #img4[j,i] += 3.0*dF*dxdy
                #img4[j,i] = 5.0

                #zipper = abs(x) < 0.18 && y > 4.67
                #if !zipper
                #spot_flux[kd] += dF
                
                #spot_flux[k] += dF
                #end
            end #if inside            

            
            end #hiti
            #end #XXX debug
        end# for x
    end#for y

    #continue
    
    #TODO: deal with hidden spot
    #i.e. skip time bin
    
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

    #frame_y1 -= frame_expansion_y2*0.15
    #frame_y2 += frame_expansion_y2*0.15
    #frame_x1 -= frame_expansion_x2*0.15
    #frame_x2 += frame_expansion_x2*0.15
    #frame_y1 *= 0.95 
    #frame_y2 *= 1.05
    #frame_x1 *= 0.95
    #frame_x2 *= 1.05

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

    #select larger
    #if frame_xs < frame_ys
    #    Nx_frame = N_frame
    #    Ny_frame = max(round(Int, (frame_ys*N_frame/frame_xs)), 2)
    #else
    #    Ny_frame = N_frame
    #    Nx_frame = max(round(Int, (frame_xs*N_frame/frame_ys)), 2)
    #end

    #keep aspect ratio
    #if frame_xs < frame_ys
    #    Nx_frame = max(round(Int, (frame_ys*N_frame/frame_xs)), 2)
    #    Ny_frame = max(round(Int, (frame_ys*N_frame/frame_xs)), 2)
    #else
    #    Ny_frame = max(round(Int, (frame_xs*N_frame/frame_ys)), 2)
    #    Nx_frame = max(round(Int, (frame_xs*N_frame/frame_ys)), 2)
    #end

    

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


    
    #Integrate flux inside of the spot image frames
    
    #img5[:,:] = 0.0
    img5 = zeros(Ny_frame, Nx_frame) #debug array

    Ndelta = 0.0

    old_subframe = [frame_y2,
                    frame_y1,
                    frame_x1,
                    frame_x2
                    ]

    for j = 1:Ny_frame
        y = frame_ygrid[j]
        for i = 1:Nx_frame
            x = frame_xgrid[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))

            if interp
                phi = phi_interp_atan(rad,chi)
                theta = theta_interp[rad,chi]
                Xob = Xs_interp[rad,chi]
                time = time_interp[rad,chi]
                cosa = cosa_interp[rad,chi]
                #hit = hits_interp[rad,chi] #test if we hit the surface

            #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
            else
                time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                            X, Osb, beta, quad, wp, Rg)
                time -= time0
            #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
            #println()
            end
            
            
            #hiti = round(Int, hit)
            #hiti = rad < rmax ? hiti : 0
            #println(hit, hiti)
            
            #if hiti > 0
            if rad <= edge_interp(chi)
                
                dt = time*G*M/c^3
                #dt = 0.0
                
                #rotatate star
                #println("t: $t dt: $dt $(dt/t)")
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
                
                img5[j,i] = painter(phi, theta)/2.0
                
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
                        old_subframe[1] = frame_y2 < y ? y : frame_y2 #top #max
                        old_subframe[2] = frame_y1 > y ? y : frame_y1 #bottom #min
                        old_subframe[3] = frame_x1 > x ? x : frame_x1 #left min
                        old_subframe[4] = frame_x2 < x ? x : frame_x2 #right max

                        
                        EEd, delta = radiation(rad, chi,
                                               phi, theta, cosa,
                                               X, Xob, Osb, sini)
                    end
                    #println(EEd, " ", delta)
                    #println()
                    
                    #println("inside")
                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)

                    #println(dfluxE)
                    #println("x=$x y=$y r=$rad chi=$chi")

                    
                    sdelta[k] += delta * frame_dxdy #* imgscale
                    sdelta2[k] += EEd * frame_dxdy #* imgscale
                    Ndelta += frame_dxdy
                    
                    img5[j,i] += dfluxB * frame_dxdy * imgscale * 1.0e5

                    #skip timelag
                    kd = k
                    
                    for ie = 1:3
                        sfluxE[kd, ie] += dfluxE[ie] * frame_dxdy * imgscale
                        sfluxNE[kd, ie] += dfluxNE[ie] * frame_dxdy * imgscale
                    end
                    sfluxNB[kd] += dfluxNB * frame_dxdy * imgscale
                    sfluxB[kd] += dfluxB * frame_dxdy * imgscale

                    #catch_NaN(sfluxE)
                    
                end #inside spot
            end#hiti
        end #x
    end#y

    println("bol flux: $(sfluxB[k]) | num flux $(sfluxNB[k])")

    p10b = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 2, "Blues")

    #add time stamp
    xs = frame_xgrid[1] + 0.84*(frame_xgrid[end]-frame_xgrid[1])
    ys = frame_ygrid[1] + 0.93*(frame_ygrid[end]-frame_ygrid[1])
    Winston.add(p10b, Winston.DataLabel(xs, ys, "$(k)"))
    #display(p10)

    #bol flux
    #p10c = plot(phase, sfluxB, "k-")
    #p10c = oplot([phase[k]], [sfluxB[k]], "ko")
    p10c = plot(phase, sfluxE[:,1], "k-")
    p10c = oplot([phase[k]], [sfluxE[k,1]], "ko")

    #doppler factor
    sdelta[k] = sdelta[k]/Ndelta
    p10d = plot(phase, sdelta, "b-")
    #p10d = oplot(phase, sdelta2, "r-")
    #p10c = oplot(phase, (sdelta2./Ndelta).^5, "r-")
    #p10c = plot(phase, (sdelta2./sdelta), "r-")
    #p10c = oplot([phase[k]], [sfluxB[k]], "ko")

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
toc()


#write to file
#opath = "out/"

#opath = "out2/"
#opath = "out2/cadeau+morsink/"
#opath = "out2/f$(round(Int,fs))/r$(round(Int,R/1e5))n/"
opath = "out3/HT/"

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
