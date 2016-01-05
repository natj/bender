#Compute spot on the NS


#Load photon bender
#include("bender.jl")

#Compute raw image
#include("comp_img.jl")

#Interpolate from raw image and compute radiation processes
#include("radiation.jl")

rho = deg2rad(1.0)
colat = deg2rad(50.0)


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
    
    return false
end


#photon time lag function
# returns new time-array index taking into account time delay effect
function time_lag(t, k, times, Nt, tbin, phi, theta)

    #exact from raytracing
    #dt = t/c
    #dt = t*G*M/c^3
    
    #approximative
    cosi = sqrt(1-sini^2)
    cospsi = cosi*cos(theta) + sini*sin(theta)*cos(phi)
    y = 1 - cospsi

    #dt0 = y*R/c
    dt = y*(1.0 + (U*y/8.0)*(1+ y*(1.0/3.0 - U/14.0)))*R/c
    #println("dt: $dt  dt0: $dt0  dt: $dt01 ")
    #println(" $(dt01/dt) ")

    #get new timebin index
    kd = 1
    while dt > (2kd - 1)*tbin
        kd += 1
    end
    kd -= 1

    kindx = k + kd
    if kindx > Nt; kindx -= Nt; end
        
    #println("k: $k kd: $kd kindx: $kindx")
    #println()
    
    return kindx
end




img4 = zeros(Ny_dense, Nx_dense) #debug array

#Spot image frame size
N_frame = 200


#Beaming function for the radiation
#Ir(cosa) = 1.0 #isotropic beaming
Ir(cosa) = cosa

#Time parameters
Nt = 64
times = collect(linspace(0, 1/fs, Nt))
tbin = abs(times[2] - times[1])/2.0 
phase = collect(times .* fs)

spot_flux = zeros(Nt)


tic()
for k = 1:Nt
#for k = 35:35
#for k = 80:80
#for k = 24:38
#for k = 20:45
        
    img4[:,:] = 0.0
    
    t = times[k]
    println("t: $t k: $k")
    
    #set whole image as starting frame
    frame_y2 = y_grid_d[1]
    frame_y1 = y_grid_d[end]
    frame_x1 = x_grid_d[end]
    frame_x2 = x_grid_d[1]

    #inf_small = true
    
    for j = y1s:y2s
        y = y_grid_d[j]

        #frame_left = false
        
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))
            
            #trace back to star
            phi = phi_interp_atan(rad,chi)
            theta = theta_interp[rad,chi]
                        
            #rotate star
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)
      
            img4[j,i] = painter(phi, theta)/2.0

            inside = spot(0.0, phi, theta,
                          stheta = colat,
                          delta = rho
                          )
            
            if inside
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
                dF = flux_interp[y,x]
                dE = reds_interp[y,x]
                                
                #img4[j,i] = painter(phi, theta)
                img4[j,i] += 3.0*dF /dxdy
                #img4[j,i] = 5.0

                #zipper = abs(x) < 0.18 && y > 4.67
                #if !zipper
                #spot_flux[kd] += dF
                
                #spot_flux[k] += dF
                #end
            end #if inside            

            

        end# for x
    end#for y

    #continue
    
    #TODO: deal with hidden spot
    #i.e. skip time bin

    
    #expand image a bit
    #########
    frame_expansion_x = abs(x_grid_d[4] - x_grid_d[1])
    frame_expansion_y = abs(y_grid_d[4] - y_grid_d[1])
    frame_y2 += frame_expansion_y
    frame_y1 -= frame_expansion_y
    frame_x1 -= frame_expansion_x
    frame_x2 += frame_expansion_x

    #Exapand image a bit keeping the aspect ratio
    ##########
    #frame_expansion_x = abs(frame_x2 - frame_x1)
    #frame_expansion_y = abs(frame_y2 - frame_y1)

    #frame_y1 -= frame_expansion_y*0.02
    #frame_y2 += frame_expansion_y*0.02
    #frame_x1 -= frame_expansion_x*0.02
    #frame_x2 += frame_expansion_x*0.02

    frame_xs = abs(frame_x2 - frame_x1)/N_frame
    frame_ys = abs(frame_y2 - frame_y1)/N_frame


    println("x1: $frame_x1  x2: $frame_x2  y1: $frame_y1 y2: $frame_y2")  
    println("x = $frame_xs y = $frame_ys")

    #pick smaller
    #if frame_xs > frame_ys
    #    Nx_frame = N_frame
    #    Ny_frame = max(round(Int, (frame_ys*N_frame/frame_xs)), 2)
    #else
    #    Ny_frame = N_frame
    #    Nx_frame = max(round(Int, (frame_xs*N_frame/frame_ys)), 2)
    #end

    if frame_xs < frame_ys
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
    p10a = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 5, "Blues")
    Winston.add(p10a, Curve([frame_x1, frame_x2, frame_x2, frame_x1, frame_x1],
                           [frame_y1, frame_y1, frame_y2, frame_y2, frame_y1],
                           linestyle="solid"))
    #add time stamp
    xs = x_grid_d[1] + 0.84*(x_grid_d[end] - x_grid_d[1])
    ys = y_grid_d[1] + 0.93*(y_grid_d[end] - y_grid_d[1])
    Winston.add(p10a, Winston.DataLabel(xs, ys, "$(k) ($(round(times[k]*fs,3)))"))
    #display(p10)

    #println("dx = $(frame_dxx) dy = $(frame_dyy)")


    
    #Integrate flux inside of the spot image frames
    
    #img5[:,:] = 0.0
    img5 = zeros(Ny_frame, Nx_frame) #debug array
    
    for j = 1:Ny_frame
        y = frame_ygrid[j]
        for i = 1:Nx_frame
            x = frame_xgrid[i]

            rad = hypot(x,y)
            chi = mod2pi(pi/2 - atan2(y,x))

            #interpolate if we are not on the edge or near the zipper
            #ring = rstar_min*0.98 < sqrt(x^2 + y^2) < 1.01*rstar_max
            #zipper = abs(x) < 0.1 && y > 3.0
            ring = false
            zipper = false
            #ring = true
            #zipper = true

            
            if ring || zipper
                time, phi, theta, Xob, hit, cosa = bender3(x, y, sini,
                                                           X, Osb,
                                                           beta, quad, wp, Rg)
            else
                # phi & theta
                phi = phi_interp_atan(rad,chi)
                theta = theta_interp[rad,chi]
                Xob = Xs_interp[rad,chi]
                time = time_interp[rad,chi]
                cosa = cosa_interp[rad,chi]
                
                #test if we hit the surface
                hit = hits_interp[rad,chi]
            end

            #time, phi, theta, Xob, hit, cosa = bender3(x, y, sini,
            #                                           X, Osb,
            #                                           beta, quad, wp, Rg)

            #test if we hit the surface
            #hit = hits_interp[y,x]
            #hiti = round(Int,hit - 0.49)
            hiti = round(Int, hit)
            
            if hiti > 0
                #phi = phi_interp_atan(y,x)
                #theta = theta_interp[y,x]

                #rotatate star
                phi = phi - t*fs*2*pi
                phi = mod2pi(phi)
                
                img5[j,i] = painter(phi, theta)/2.0
                #if (ring || zipper)
                #    img5[j,i] = painter(phi, theta)/2.0
                #end

                
                inside = spot(0.0, phi, theta,
                              stheta = colat,
                              delta = rho
                              )

                if inside
                    #time = time_interp[y,x]

                    #compute 
                    #earea = polyarea(x, y,
                    #                 frame_dxx, frame_dyy,
                    #                 phi, theta,
                    #                 exact=(ring || zipper)
                    #                 #exact=true
                    #                 #exact=false
                    #                 )

                    #kd = time_lag(time, k, times, Nt, tbin, phi, theta)

                    #Xob = Xs_interp[y,x] 
                    #cosa = cosa_interp[y,x]
                    #dF, dE = radiation(Ir,
                    #                   x,y,
                    #                   phi, theta, cosa,
                    #                   X, Xob, Osb, sini, earea)
                    
                    dF = flux_interp[y,x]
                    dE = reds_interp[y,x]
                    
                    
                    #img4[j,i] = painter(phi, theta)
                    
                    #img5[j,i] += 1.0e9*dF * frame_dxdy
                    #println(dF)
                    img5[j,i] += 1.0e5 * dF
                    
                    #img5[j,i] = 5.0
                    #spot_flux[kd] += dF * frame_dxdy
                    #spot_flux[kd] += frame_dxdy
                    spot_flux[k] += dF * frame_dxdy
                end #inside spot

                
            end#hiti
        end #x
    end#y

    p10b = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 5, "Blues")

    #add time stamp
    xs = frame_xgrid[1] + 0.84*(frame_xgrid[end]-frame_xgrid[1])
    ys = frame_ygrid[1] + 0.93*(frame_ygrid[end]-frame_ygrid[1])
    Winston.add(p10b, Winston.DataLabel(xs, ys, "$(k)"))
    #display(p10)

    p10c = plot(phase, spot_flux, "k-")
    p10c = oplot([phase[k]], [spot_flux[k]], "ko")

    tt = Table(3,1)
    tt[1,1] = p10a
    tt[2,1] = p10b
    tt[3,1] = p10c
    display(tt)

    #readline(STDIN)

end#for t
toc()


#write to file
opath = "out/"
mkpath(opath)

#fname = "f$(fs)_lamb_bb_R$(round(R/1e5,1))_M$(round(M/Msun,1))_rho30.csv"
fname = "f$(round(Int,fs))_lamb_bb_R$(round(R/1e5,1))_M$(round(M/Msun,1))_rho$(round(Int,rad2deg(rho))).csv"
wmatr = zeros(Nt, 2)
wmatr[:,1] = phase
wmatr[:,2] = spot_flux
writecsv(opath*fname, wmatr)
