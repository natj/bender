#Compute spot on the NS

rho = deg2rad(30.0)
colat = deg2rad(50.0)

#rho = deg2rad(10.0)
#colat = deg2rad(90.0)

interp = true
exact_edges = false


function spotr(t, phi, theta,
                       stheta, delta,
                       X, Xob, Osb)

    nu2   = beta/3.0 - quad*0.5*(3*cos(theta)^2-1)
    B2    = beta
    zeta2 = beta*((4/3)*0.5*(3*cos(theta)^2 - 1) - 1/3)
    Rgm, dR = Rgmf(theta, X, Osb)

    enu = (1-Xob/2)/(1+Xob/2)*exp(nu2*Xob^3)
    B = (1-Xob/2)*(1+Xob/2) + B2*Xob^2
    ezeta = (1-Xob/2)*(1+Xob/2)*exp(zeta2*Xob^2)

    w = wp*Xob^3*(1-3*Xob) /(G*M/c^3) #into rad/seconds
    vz = Rgm*(1/enu)*sin(theta)*(2pi*fs - w) #isoradial zamo
    bz = R*vz/c
    gamma = 1/sqrt(1 - bz^2)
    #gamma = 1.0

    #if (t/gamma) != 0
    #    println("gamma comparison $(t/gamma)")
    #end

    const rho = 1.0
    bet=ezeta/B
    the = theta

    #spot center
    phi1 = 0.0
    the1 = stheta
    x1 = cos(phi1)*sin(the1)
    y1 = sin(phi1)*sin(the1)
    z1 = cos(the1)
    

    #photon location
    phi2 = phi
    the2 = theta
    x2 = cos(phi2)*sin(the2)
    y2 = sin(phi2)*sin(the2)
    z2 = cos(the2)
    
    #dot product (O(Omega^2) in rotation, i.e. only y-dir contracted)
    r1r2= x1*x2 + gamma*gamma*y1*y2 + z1*z2
    r1r1 = x1^2 + gamma*gamma*y1^2 + z1^2 
    r2r2 = x2^2 + gamma*gamma*y2^2 + z2^2 

    d = acos(r1r2/sqrt(r1r1)/sqrt(r2r2))

    if abs(d) <= delta
        return true
    end

    return false
end


function spot(gamma, phi, theta;
              stheta = deg2rad(50.0), #spot latitude
              delta = deg2rad(30.0) #spot half-opening angle
              )

    #Circular spot
    d = great_circle_dist_rel(0.0, phi, stheta, theta, gamma)
    #d = great_circle_dist(0.0, phi, stheta, theta)
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
N_frame = 50

N_frame_chi = 500
N_frame_rad = 100

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
sdelta3 = zeros(Nt)
itersv = zeros(Nt)

sfluxE = zeros(Nt, 3)
sfluxNE = zeros(Nt, 3)

sfluxB = zeros(Nt)
sfluxNB = zeros(Nt)

iters = zeros(1)

function dFcart(xx, yy)
    x = xx[1]
    y = xx[2]
    ##println("dF: x:$x y:$y $imgscale $t")

    rad = hypot(x,y)
    chi = mod2pi(pi/2 - atan2(y,x))
    
    dS = 1.0

    return dF(rad, chi, dS, yy)
end

function dFpol(xx, yy)
    rad = xx[1]
    chi = mod2pi(xx[2])

    #println("rad $rad chi $chi")

    dS = rad
    
    return dF(rad, chi, dS, yy)
end


function dF(rad, chi, dS, yy)

    dfluxNB = 0.0
    dfluxB = 0.0
    dfluxE = zeros(3)
    dfluxNE = zeros(3)

    if interp
        phi = phi_interp_atan(rad,chi)
        theta = theta_interp[rad,chi]
        Xob = Xs_interp[rad,chi]
        time = time_interp[rad,chi]
        cosa = cosa_interp[rad,chi]
        #dtau = dtau_interp[rad,chi]
        #hit = hits_interp[rad,chi] #test if we hit the surface
    
        #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
    else
        time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                    X, Osb, beta, quad, wp, Rg)
        time -= time0
        #println(phi," ",theta," ",Xob," ",time," ",cosa," ",hit)
        #println()
    end

    if rad <= edge_interp(chi) -1.0e-2
        dt = time*G*M/c^3
        phi = phi - (t - dt)*fs*2*pi
        phi = mod2pi(phi)

        EEd, delta, gamma = radiation(rad, chi,
                                     phi, theta, cosa,
                                     X, Xob, Osb, sini)

        inside = spotr(0.0, phi, theta,
                       colat, rho,
                       X, Xob, Osb)
        #inside = spot(gamma, phi, theta,
        #              stheta = colat,
        #              delta = rho
        #             )
        if inside
    
            if interp
                unity = unity_interp[rad,chi]
                #EEd = reds_interp[rad,chi]
                #println(EEd, " ", delta)
            else
    
                # update subimage corners; they might no be up-to-date
                # because we trace exactly here instead of interpolating.
                #old_subframe[1] = frame_y2 < y ? y : frame_y2 #top #max
                #old_subframe[2] = frame_y1 > y ? y : frame_y1 #bottom #min
                #old_subframe[3] = frame_x1 > x ? x : frame_x1 #left min
                #old_subframe[4] = frame_x2 < x ? x : frame_x2 #right max
    
                
                #EEd, delta, dtau = radiation(rad, chi,
                #                             phi, theta, cosa,
                #                             X, Xob, Osb, sini)
            end
            #println(EEd, " ", delta)
            #println()
   	    
            if (exact_edges && unity != 1.0) || (exact_edges && rad <= edge_interp(chi) - 1.0e-3)
                #println("exact edge for $rad")
                #print(".")
                time, phi, theta, Xob, hit, cosa = bender3p(rad, chi, sini,
                                                            X, Osb, beta, quad, wp, Rg)
                time -= time0
                dt = time*G*M/c^3
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
                
		EEd, delta, gamma = radiation(rad, chi,
                                                 phi, theta, cosa,
                                                 X, Xob, Osb, sini)
                #inside = spot(gamma, phi, theta,
                #              stheta = colat,
                #              delta = rho)
                inside = spotr(0.0, phi, theta,
                       colat, rho,
                       X, Xob, Osb)

                if inside && hit
   			#exact edge 
                else
                    EEd = 0.0
                    delta = 0.0
                    dtau = 0.0
                    gamma = 1.0
                end
            end
    
            dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)

    
            #dfluxB *= 1/gamma
            #dfluxNB *= 1/gamma
            #for ie = 1:3
            #   dfluxE[ie] *= 1/gamma 
            #   dfluxNE[ie] *= 1/gamma
            #end
    
        end #inside spot
    end#hiti

    #println(dfluxB)
    yy[1] = dS*dfluxE[1]
    yy[2] = dS*dfluxE[2]
    yy[3] = dS*dfluxE[3]

    yy[4] = dS*dfluxNE[1]
    yy[5] = dS*dfluxNE[2]
    yy[6] = dS*dfluxNE[3]

    yy[7] = dS*dfluxNB
    yy[8] = dS*dfluxB

    #yy[:] = [dfluxE[1],dfluxE[2],dfluxE[3],
    #         dfluxNE[1],dfluxNE[2],dfluxNE[3],
    #         dfluxNB,
    #         dfluxB
    #        ]

    iters[1] += 1

    return 1
end


###################################
#Polar integration

#create thick radius grid
#Nrad_frame = 1000
#rad_diffs = 1 ./ exp(linspace(0.0, 1.5, Nrad_frame-1).^2)
#rad_grid_d = rmax * cumsum(rad_diffs) / sum(rad_diffs)
#unshift!(rad_grid_d, 0.0)

chi_sweep = zeros(Nchi)

old_subframe = [y_grid_d[1],
                y_grid_d[end],
                x_grid_d[end],
                x_grid_d[1]
                ]

t = 0.0
frame_dxdy = 1.0

tic()
for k = 1:Nt
#for k = 45:65
#for k = 2:1
    t = times[k]
    println(" ")
    println("t: $t k: $k phase $(phase[k])")

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
                Xob = Xs_interp[rad,chi]

                #rotate star
                dt = time*G*M/c^3
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)

                inside = spotr(0.0, phi, theta,
                       colat, rho,
                       X, Xob, Osb)

                #inside = spot(0.0, phi, theta,
                #              stheta = colat,
                #              delta = rho
                #              )

            
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

    if !located_spot
    	continue
    end
    
    #bounding box limits
    jradmin = max(1, jradmin - 4)
    jradmax = min(Nrad, jradmax + 4)
    println("radmin: $(rad_grid[jradmin]) radmax: $(rad_grid[jradmax])")
    #println("chimin: $(chi_grid[ichimin]) chimax: $(chi_grid[ichimax])")

    #build chi limits for bounding box
    polar_integration = true
    if jradmin == 1 #full circle
        polar_integration = false

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
    
    #make pizza slize bigger
    #chimin -= 0.01
    #chimax += 0.01

    
    #plot in cartesian coords
    img4[:,:] = 0.0

    frame_y2 = y_grid_d[1]
    frame_y1 = y_grid_d[end]
    frame_x1 = x_grid_d[end]
    frame_x2 = x_grid_d[1]
    
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
                cosa = cosa_interp[rad,chi]
                Xob = Xs_interp[rad,chi]
            
                #rotate star
                dt = time*G*M/c^3 #time shift
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
      
        	EEd, delta, gamma = radiation(rad, chi,
                                     phi, theta, cosa,
                                     X, Xob, Osb, sini)

                img4[j,i] = painter(phi, theta)/2.0
                
                inside = spotr(0.0, phi, theta,
                       colat, rho,
                       X, Xob, Osb)

                #inside = spot(0.0, phi, theta,
                #              stheta = colat,
                #              delta = rho
                #              )
             
                if inside
                    frame_y2 = frame_y2 < y ? y : frame_y2 #top #max
                    frame_y1 = frame_y1 > y ? y : frame_y1 #bottom #min
                    frame_x1 = frame_x1 > x ? x : frame_x1 #left min
                    frame_x2 = frame_x2 < x ? x : frame_x2 #right max

                    dfluxE, dfluxNE, dfluxNB, dfluxB = bbfluxes(EEd, delta, cosa)
                    img4[j,i] += dfluxB * dxdy * imgscale*1.0e7
                end #if inside            
            end #hiti
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


    #Cartesian subgrid integration
    img5 = zeros(Ny_frame, Nx_frame) #subgrid img
    img5[1,1] = 1.0e-5 #prevent plotting from crashing with empty array
    Ndelta = 0.0

    old_subframe = [frame_y2,
                    frame_y1,
                    frame_x1,
                    frame_x2
                    ]

    iters[1] = 0
    for j = 1:Ny_frame
        y = frame_ygrid[j]
        for i = 1:Nx_frame
            x = frame_xgrid[i]

            dFs = zeros(8)
            istat = dFcart([x,y], dFs) 
            dFs *= imgscale * frame_dxdy

            sfluxE[k, 1] += dFs[1]
            sfluxE[k, 2] += dFs[2]
            sfluxE[k, 3] += dFs[3]

            sfluxNE[k, 1] += dFs[4]
            sfluxNE[k, 2] += dFs[5]
            sfluxNE[k, 3] += dFs[6]

            sfluxNB[k] += dFs[7]
            sfluxB[k] += dFs[8]

            img5[j,i] += dFs[8] * frame_dxdy * imgscale * 5.0e5
        end #x
    end#y
    println("Riemann integral iterations: $(iters[1])")
    println("num flux $(sfluxNB[k])")

    #plot subcartesian grid
    p10b = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 0, "Blues")

    #add time stamp
    xs = frame_xgrid[1] + 0.84*(frame_xgrid[end]-frame_xgrid[1])
    ys = frame_ygrid[1] + 0.93*(frame_ygrid[end]-frame_ygrid[1])
    Winston.add(p10b, Winston.DataLabel(xs, ys, "$(k)"))
    #display(p10)

    Winston.add(p10b, Curve(frame_xs, frame_ys,
                            linestyle="solid"))
    Winston.add(p10b, Curve([frame_xs[1], frame_xs2[1]],
                            [frame_ys[1], frame_ys2[1]],
                            linestyle="solid",
                            color="black"))
    Winston.add(p10b, Curve([frame_xs[end], frame_xs2[end]],
                            [frame_ys[end], frame_ys2[end]],
                            linestyle="solid",
                            color="red"))
    Winston.add(p10b, Curve(frame_xs2, frame_ys2,
                            linestyle="solid", color="red"))




    #linear chi grid with thick spacing
    #chi_grid2 = linspace(chimin, chimax, N_frame_chi)
    #chi_diff2 = zeros(N_frame_chi)
    #
    #chi_diff2[2:N_frame_chi] = 0.5 * abs(diff(chi_grid2))
    #chi_diff2[1:N_frame_chi-1] .+= 0.5 * abs(diff(chi_grid2))

    #chi_diff2[1] = 0.5 * abs(chi_grid2[2] - chi_grid[1])
    #chi_diff2[N_frame_chi] = 0.5 * abs(chi_grid2[N_frame_chi] - chi_grid2[N_frame_chi-1])

    ##weighted rad grid with double spacing
    #rad_diff2 = zeros(Nrad)

    #rad_diff2[2:Nrad] = 0.5 * abs(diff(rad_grid))
    #rad_diff2[1:Nrad-1] .+= 0.5 * abs(diff(rad_grid))

    #rad_diff2[1] = 0.5 * abs(rad_grid[2] - rad_grid[1])
    #rad_diff2[Nrad] = 0.5 * abs(rad_grid[Nrad] - rad_grid[Nrad-1])
    #println(chi_grid2)

    #for i = 1:N_frame_chi
    #    chi = mod2pi(chi_grid2[i])
    #    for j = jradmin:jradmax
    #        rad = rad_grid[j]
    #        end
    #    end
    #end



    tic()
    iters[1] = 0
    vals = zeros(8)
 
    if located_spot 
    if polar_integration
    #polar cubature integration
    println(" Polar Cubature integration")

    (vals, errs) = pcubature(8,
    #(vals, errs) = hcubature(8,
                            dFpol,
                            [rad_grid[jradmin], chimin],
                            [rad_grid[jradmax], chimax],
                            reltol = 1.0e-4,
                            abstol = 0.0,
                            #maxevals=int(1e5))
                            maxevals=int(2e7))
                            #maxevals=int(2e6))

    else #cartesian
    println(" Cartesian Cubature integration")
    (vals, errs) = pcubature(8,
                           dFcart,
                           [frame_x1, frame_y1],
                           [frame_x2, frame_y2],
                           reltol = 1.0e-4,
                           abstol = 0.0,
                           maxevals=int(2e7))

    end #if-else polar/cartesian
    println("Cubature iterations: $(iters[1])")
    toc()
    end

    vals *= imgscale
    errs = errs ./ vals #transform into relative err
    errs *= imgscale

    val = vals[7]
    err = errs[7]
    #merr = (sfluxNB[k] - val)/val

    sfluxE[k, 1] = vals[1]
    sfluxE[k, 2] = vals[2]
    sfluxE[k, 3] = vals[3]

    sfluxNE[k, 1] = vals[4]
    sfluxNE[k, 2] = vals[5]
    sfluxNE[k, 3] = vals[6]

    sfluxNB[k] = vals[7]
    sfluxB[k] = vals[8]

    println("num flux $val $err")


    #plot
    #bol flux
    p10c = plot(phase, sfluxB, "k-")
    p10c = oplot([phase[k]], [sfluxB[k]], "ko")

    #doppler factor
    sdelta[k] = sdelta3[k]/Ndelta
    sdelta2[k] = sdelta3[k]/Ndelta
    sdelta3[k] = sdelta3[k]/Ndelta
    sdelta3[k] = Ndelta

    #p10d = plot(phase, sdelta3, "b-")
    itersv[k] = iters[1]
    p10d = plot(phase, itersv, "b-")

    tt1 = Table(1,2)
    tt1[1,1] = p10a
    tt1[1,2] = p10b

    tt = Table(3,1)
    #tt[1,1] = p10a
    tt[1,1] = tt1
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

t = 0.0
frame_dxdy = 1.0


#for k = 1:Nt
for k = 2:1
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
                cosa = cosa_interp[rad,chi]
        	Xob = Xs_interp[rad,chi]
                
                #rotate star
                dt = time*G*M/c^3 #time shift
                #dt = 0.0
                
                phi = phi - (t - dt)*fs*2*pi
                phi = mod2pi(phi)
      
                #img4[j,i] = -4.0*painter(phi, theta)/2.0
                img4[j,i] = painter(phi, theta)/2.0
                
        	EEd, delta, gamma = radiation(rad, chi,
                                     phi, theta, cosa,
                                     X, Xob, Osb, sini)

                inside = spotr(0.0, phi, theta,
                       colat, rho,
                       X, Xob, Osb)
                #inside = spot(gamma, phi, theta,
                #              stheta = colat,
                #              delta = rho
                #              )
             
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
                    #time = time_interp[rad,chi]
                    #cosa = cosa_interp[rad,chi]
                    #delta = delta_interp[rad,chi]
                    #EEd = reds_interp[rad,chi]
                        
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



    img5 = zeros(Ny_frame, Nx_frame) #subgrid img

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
            istat = dFcart([x,y], dFs) 
            dFs *= imgscale * frame_dxdy

            sfluxE[k, 1] += dFs[1]
            sfluxE[k, 2] += dFs[2]
            sfluxE[k, 3] += dFs[3]

            sfluxNE[k, 1] += dFs[4]
            sfluxNE[k, 2] += dFs[5]
            sfluxNE[k, 3] += dFs[6]

            sfluxNB[k] += dFs[7]
            sfluxB[k] += dFs[8]

            img5[j,i] += dFs[8] * frame_dxdy * imgscale * 5.0e5
        end #x
    end#y
    println("Riemann integral iterations: $(iters[1])")
    toc()

    #println("bol flux: $(sfluxB[k]) | num flux $(sfluxNB[k])")
    println("num flux $(sfluxNB[k])")

    #start cubature integration only if we see the spot

    val = 0.0
    if (sfluxNB[k] > 0)

    tic()
    iters[1] = 0
    (vals, errs) = pcubature(8,
                           dFcart,
                           [frame_x1, frame_y1],
                           [frame_x2, frame_y2],
                           reltol = 1.0e-4,
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

    end #end of cubature integration

    println("num flux $val $err $merr")
    p10b = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 0, "Blues")

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

    #p10d = plot(phase, sdelta, "b-",
    #            xrange=[0.0, 1.0],
    #            yrange=[0.8, 1.2])

    itersv[k] = iters[1]
    p10d = plot(phase, itersv, "b-")



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
#opath = "out/"
#opath = "out2/"
#opath = "out2/cadeau+morsink/"
#opath = "out2/f$(round(Int,fs))/r$(round(Int,R/1e5))n/"
opath = "out3/"
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



