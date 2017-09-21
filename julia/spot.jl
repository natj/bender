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
exact_edges = true


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
N_frame = 100

N_frame_chi = 500
N_frame_rad = 100

#Beaming function for the radiation
#Ir(cosa) = 1.0 #isotropic beaming
#Ir(cosa) = cosa

#Time parameters
Nt = 32

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

#Polar integration

#create thick radius grid
Nrad_frame = 1000
rad_diffs = 1 ./ exp(linspace(0.0, 1.5, Nrad_frame-1).^2)
rad_grid_d = rmax * cumsum(rad_diffs) / sum(rad_diffs)
unshift!(rad_grid_d, 0.0)

chi_sweep = zeros(Nchi)

tic()
#for k = 1:Nt
for k = 2:1
    t = times[k]

    jradmin = Nrad
    jradmax = 1
    ichimin = Nchi
    ichimax = 1

    located_spot = false
    
    for i = 2:Nchi-1
        chi = chi_grid[i]

        #avoid double integrating angle edges
        #if chi <= 0; continue; end
        #if chi >= 2pi; continue; end
        
        #if mod(i,10) == 0; println("chi: $(round(chi/2pi,2))"); end

        located_spot_chi = false
        for j = 2:Nrad-1
            rad = rad_grid[j]

            if rad <= edge_interp(chi)
                
                #Ray traced photons
                phi = Phis[j, i]
                theta = Thetas[j, i]
                time = Times[j, i]
                
                #rotate star
                dt = time*G*M/c^3
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)

                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho
                              )

            
                if inside
                    located_spot = true
                    located_spot_chi = true
                    
                    #println("inside")
                    jradmax = rad > rad_grid[jradmax] ? j : jradmax
                    jradmin = rad < rad_grid[jradmin] ? j : jradmin
                    
                    #time = Times[j, i]
                    #Xob = Xs[j, i]
                    #cosa = cosas[j, i]
                    #dF = dFlux[j, i]
                    #dE = Reds[j, i]


                    #drdchi = abs(rad_grid[j+1] - rad_grid[j-1])*abs(chi_grid[i+1] - chi_grid[i-1])*rad_grid[j]
                    #spot_flux2[kd] += dF * drdchi
                end#inside spot
            end#inside star
        end#rad

        chi_sweep[i] = located_spot_chi ? 1.0 : 0.0
    end#chi

    
    #bounding box limits
    #jradmin = jradmin == 1 ? 1 : jradmin - 1
    #jradmax = jradmax == Nrad ? Nrad : jradmax + 1
    jradmin = max(1, jradmin - 2)
    jradmax = min(Nrad, jradmax + 2)
    println("radmin: $(rad_grid[jradmin]) radmax: $(rad_grid[jradmax])")
    #println("chimin: $(chi_grid[ichimin]) chimax: $(chi_grid[ichimax])")

    #build chi limits for bounding box
    if jradmin == 1 #full circle
        #ichimax = Nchi-1
        #ichimin = 2
        chimin = 0.0
        chimax = 2.0pi
    else #pizza slices
        #println(chi_sweep)

        if chi_sweep[Nchi-1] == 1 #splitted spot 
            #ichimin = minimum(findin(chi_sweep, [1]))
            #ichimax = maximum(findin(chi_sweep, [1]))

            #get end of the spot sweep
            q1 = 2
            while chi_sweep[q1] == 1
                q1 += 1
            end
            q1 += 1
            
            #get start of the spot sweep
            q2 = Nchi-1
            while chi_sweep[q2] == 1
                q2 -= 1
            end
            q2 -= 1

            chimin = chi_grid[q2] - 2.0pi
            chimax = chi_grid[q1]
            
        else
            #find limits
            ichimin = minimum(findin(chi_sweep, [1]))
            ichimax = maximum(findin(chi_sweep, [1]))

            #expand
            ichimin = ichimin == 1 ? 1 : ichimin - 1
            ichimax = ichimax == Nchi ? Nchi : ichimax + 1

            chimin = chi_grid[ichimin]
            chimax = chi_grid[ichimax]
        end
    end
    

    
    #plot in cartesian coords
    img4[:,:] = 0.0
    
    for j = y1s:y2s
        y = y_grid_d[j]
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))

            if rad <= edge_interp(chi)
            
                #trace back to star
                phi = phi_interp_atan(rad,chi)
                theta = theta_interp[rad,chi]
                time = time_interp[rad,chi]
            
                #rotate star
                dt = time*G*M/c^3 #time shift
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
      
                img4[j,i] = painter(phi, theta)/2.0
                
                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho
                              )
             
                if inside #&& y > -5

                    cosa = cosa_interp[rad,chi]
                    delta = delta_interp[rad,chi]
                    EEd = reds_interp[rad,chi]
                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)
                    img4[j,i] += dfluxB * dxdy * imgscale*1.0e7
                end #if inside            
            end #hiti
        end# for x
    end#for y
    
    
    #Plot large image with bounding box for the spot
    p10a = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 2, "Blues")
    #add time stamp
    xs = x_grid_d[1] + 0.84*(x_grid_d[end] - x_grid_d[1])
    ys = y_grid_d[1] + 0.93*(y_grid_d[end] - y_grid_d[1])
    Winston.add(p10a, Winston.DataLabel(xs, ys, "$(k) ($(round(times[k]*fs,3)))"))

    frame_Nchi = 20
    frame_chis = linspace(chimin, chimax, frame_Nchi)
    frame_xs = zeros(frame_Nchi)
    frame_ys = zeros(frame_Nchi)
    frame_xs2 = zeros(frame_Nchi)
    frame_ys2 = zeros(frame_Nchi)
    
    for i = 1:frame_Nchi
        chi = frame_chis[i]

        #transform to cartesian
        frame_xs[i] = rad_grid[jradmax]*sin(chi) 
        frame_ys[i] = rad_grid[jradmax]*cos(chi)

        frame_xs2[i] = rad_grid[jradmin]*sin(chi) 
        frame_ys2[i] = rad_grid[jradmin]*cos(chi)
    end
        
    Winston.add(p10a, Curve(frame_xs, frame_ys,
                            linestyle="solid"))
    Winston.add(p10a, Curve([frame_xs[1], frame_xs2[1]],
                            [frame_ys[1], frame_ys2[1]],
                            linestyle="solid",
                            color="black"))
    Winston.add(p10a, Curve([frame_xs[end], frame_xs2[end]],
                            [frame_ys[end], frame_ys2[end]],
                            linestyle="solid",
                            color="red"))
    Winston.add(p10a, Curve(frame_xs2, frame_ys2,
                            linestyle="solid", color="red"))
    
    #integrate
    Ndelta = 0.0

    #linear chi grid with thick spacing
    chi_grid2 = linspace(chimin, chimax, N_frame_chi)
    chi_diff2 = zeros(N_frame_chi)
    
    chi_diff2[2:N_frame_chi] = 0.5 * abs(diff(chi_grid2))
    chi_diff2[1:N_frame_chi-1] .+= 0.5 * abs(diff(chi_grid2))

    chi_diff2[1] = 0.5 * abs(chi_grid2[2] - chi_grid[1])
    chi_diff2[N_frame_chi] = 0.5 * abs(chi_grid2[N_frame_chi] - chi_grid2[N_frame_chi-1])

    #weighted rad grid with double spacing
    rad_diff2 = zeros(Nrad)

    rad_diff2[2:Nrad] = 0.5 * abs(diff(rad_grid))
    rad_diff2[1:Nrad-1] .+= 0.5 * abs(diff(rad_grid))

    rad_diff2[1] = 0.5 * abs(rad_grid[2] - rad_grid[1])
    rad_diff2[Nrad] = 0.5 * abs(rad_grid[Nrad] - rad_grid[Nrad-1])
    println(chi_grid2)

    for i = 1:N_frame_chi
        chi = mod2pi(chi_grid2[i])

        for j = jradmin:jradmax
            rad = rad_grid[j]

            if rad <= edge_interp(chi)

                #trace back to star
                phi = phi_interp_atan(rad,chi)
                theta = theta_interp[rad,chi]
                time = time_interp[rad,chi]
                
                

                #rotate star
                dt = time*G*M/c^3 #time shift
                #dt = 0.0
            
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)

                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho
                              )

                #if inside spot
                if inside
                    cosa = cosa_interp[rad,chi]
                    delta = clamp(delta_interp[rad,chi], 0.0, Inf)
                    EEd = clamp(reds_interp[rad,chi], 0.0, Inf)
                    dtau = clamp(dtau_interp[rad,chi], 0.0, Inf)
                    

                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)

                    #drdchi = rad * (rad_diff2[j-1] + rad_diff2[j]) * (chi_diff2[i-1] + chi_diff2[i])
                    #drdchi = 0.25*abs(rad_grid[j+1] - rad_grid[j-1])*abs(chi_grid2[i+1] - chi_grid2[i-1])*rad
                    drdchi = rad * rad_diff2[j] * chi_diff2[i]
                                        
                    sdelta[k] += delta * drdchi #* imgscale
                    sdelta2[k] += EEd * drdchi #* imgscale
                    sdelta3[k] += dtau * drdchi
                    Ndelta += drdchi
                          
                    for ie = 1:3
                        sfluxE[k, ie] += dfluxE[ie] * drdchi * imgscale *dtau
                        sfluxNE[k, ie] += dfluxNE[ie] * drdchi * imgscale *dtau
                    end
                    sfluxNB[k] += dfluxNB * drdchi * imgscale *dtau
                    sfluxB[k] += dfluxB * drdchi * imgscale *dtau

                    Ndelta += drdchi

                    #if sfluxB[k] < 0
                    #    println(cosa)
                    #    println(delta)
                    #    println(EEd)
                    #    println(dtau)
                    #end
                    
                end

            end
        end
    end

    #plot
    #bol flux
    p10c = plot(phase, sfluxB, "k-")
    p10c = oplot([phase[k]], [sfluxB[k]], "ko")

    #doppler factor
    sdelta[k] = sdelta3[k]/Ndelta
    sdelta2[k] = sdelta3[k]/Ndelta
    sdelta3[k] = sdelta3[k]/Ndelta
    sdelta3[k] = Ndelta

    p10d = plot(phase, sdelta3, "b-")

    tt = Table(3,1)
    tt[1,1] = p10a
    tt[2,1] = p10c
    tt[3,1] = p10d
    display(tt)
    
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

#for k = 2:1
for k = 1:Nt
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
             
            if inside #&& y > -5
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
                #kd = time_lag(time, k, times, Nt, tbin, phi, theta)
    
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
            
            
            #hiti = round(Int, hit)
            #hiti = rad < rmax ? hiti : 0
            #println(hit, hiti)
            
            #if hiti > 0
            if rad <= edge_interp(chi)

                #time delay
                #cosi = sqrt(1-sini^2)
                #cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
                #yparam = 1 - cospsi
                #dt = yparam*(1.0 + (U*yparam/8.0)*(1+ yparam*(1.0/3.0 - U/14.0)))*R/c
                
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

                if inside #&& y > -5

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

                        
                        EEd, delta, dtau = radiation(rad, chi,
                                                     phi, theta, cosa,
                                                     X, Xob, Osb, sini)
                    end
                    #println(EEd, " ", delta)
                    #println()

                    if exact_edges && delta < 0.98
                        println("exact edge for $rad")
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
                    
                    #println("inside")
                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)

                    #println(dfluxE)
                    #println("x=$x y=$y r=$rad chi=$chi")

                    
                    sdelta[k] += delta * frame_dxdy #* imgscale
                    sdelta2[k] += EEd * frame_dxdy #* imgscale
                    sdelta3[k] += dtau * frame_dxdy
                    Ndelta += frame_dxdy
                    
                    img5[j,i] += dfluxB * frame_dxdy * imgscale * 1.0e5

                    #skip timelag
                    #kd = k
                    
                    for ie = 1:3
                        sfluxE[k, ie] += dfluxE[ie] * frame_dxdy * imgscale *dtau
                        sfluxNE[k, ie] += dfluxNE[ie] * frame_dxdy * imgscale *dtau
                    end
                    sfluxNB[k] += dfluxNB * frame_dxdy * imgscale *dtau
                    sfluxB[k] += dfluxB * frame_dxdy * imgscale *dtau
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
    p10c = plot(phase, sfluxNB, "k-")
    p10c = oplot([phase[k]], [sfluxNB[k]], "ko")
    #p10c = plot(phase, sfluxE[:,1], "k-")
    #p10c = oplot([phase[k]], [sfluxE[k,1]], "ko")

    #doppler factor
    sdelta[k] = sdelta3[k]/Ndelta
    sdelta2[k] = sdelta3[k]/Ndelta
    sdelta3[k] = sdelta3[k]/Ndelta

    p10d = plot(phase, sdelta, "b-",
                xrange=[0.0, 1.0],
                yrange=[0.8, 1.2])
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
opath = "out/"
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
