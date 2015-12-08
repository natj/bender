#Compute spot on the NS


#Load photon bender
#include("bender.jl")

#Compute raw image
#include("comp_img.jl")

#Interpolate from raw image
#include("img.jl")



########
function spot(t, phi, theta;
              stheta = deg2rad(50.0), #spot latitude
              delta = 1.0*pi/180 #spot half-opening angle
              )

    #Vincenty's formula
    d = great_circle_dist(0.0, phi, stheta, theta)
    
    if abs(d) < delta
        return true
    end
    
    return false
end

img4 = zeros(Ny_dense, Nx_dense) #debug array

#Spot image frame size
Nx_frame = 100
Ny_frame = 100
img5 = zeros(Ny_frame, Nx_frame) #debug array



#Beaming function for the radiation
Ir(cosa) = 1.0 #isotropic beaming


#Time parameters
Nt = 64
times = linspace(0, 1/fs, Nt)
spot_flux = zeros(Nt)


tic()
for k = 1:Nt
#for k = 90:90
#for k = 80:80
#for k = 24:38
#for k = 20:45
        
    img4[:,:] = 0.0
    
    t = times[k]

    
    #set whole image as starting frame
    frame_y2 = y_grid_d[1]
    frame_y1 = y_grid_d[end]
    frame_x1 = x_grid_d[end]
    frame_x2 = x_grid_d[1]

    
    for j = y1s:y2s
        y = y_grid_d[j]

        #frame_left = false
        
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            theta = theta_interp[y,x]

            #rotate star
            phi = phi_interp_atan(y,x)
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)
      
            img4[j,i] = painter(phi, theta)/2.0

            inside = spot(0.0, phi, theta)
            
            if inside
                
                #track down spot edges
                frame_y2 = frame_y2 < y ? y : frame_y2 #top #max
                frame_y1 = frame_y1 > y ? y : frame_y1 #bottom #min
                frame_x1 = frame_x1 > x ? x : frame_x1 #left min
                frame_x2 = frame_x2 < x ? x : frame_x2 #right max
                
                
                time = time_interp[y,x]
                
                #Xob = Xs_interp[y,x] 
                #cosa = cosa_interp[y,x]
                #dF, dE = radiation(Ir,
                #                   x,y,
                #                   phi, theta, cosa,
                #                   X, Xob, Osb, sini, img3[j,i])
                dF = flux_interp[y,x]
                dE = reds_interp[y,x]
                                
                #img4[j,i] = painter(phi, theta)
                img4[j,i] += 3.0*dF/dxdy

                #img4[j,i] = 5.0
                #spot_flux[k] += dF
            end #if inside            

            

        end# for x
    end#for y

    #TODO: deal with hidden spot
    
    #expand image a bit
    frame_expansion = abs(x_grid_d[2] - x_grid_d[1])*5.0
    frame_y2 += frame_expansion
    frame_y1 -= frame_expansion
    frame_x1 -= frame_expansion
    frame_x2 += frame_expansion
    
    p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 5, "Blues")
    Winston.add(p10, Curve([frame_x1, frame_x2, frame_x2, frame_x1, frame_x1],
                           [frame_y1, frame_y1, frame_y2, frame_y2, frame_y1],
                           linestyle="dashed"))
    
    #display(p10)

    #Integrate flux inside of the spot image frames
    frame_xgrid = collect(linspace(frame_x1, frame_x2, Nx_frame))
    frame_ygrid = collect(linspace(frame_y1, frame_y2, Nx_frame))
    frame_dxx = 2.*(frame_xgrid[2] - frame_xgrid[1])
    frame_dyy = 2.*(frame_ygrid[2] - frame_ygrid[1])
    frame_dxdy = frame_dxx*frame_dyy #*X^2

    println("dx = $(frame_dxx) dy = $(frame_dyy)")
    
    img5[:,:] = 0.0
    
    for j = 1:Ny_frame
        y = frame_ygrid[j]
        for i = 1:Nx_frame
            x = frame_xgrid[i]

            
            #time, phi, theta, Xob, hit, cosa = bender3(x, y, sini,
            #                                           X, Osb,
            #                                           beta, quad, wp, Rg)

            #test if we hit the surface
            hit = hits_interp[y,x]
            hiti = round(Int,hit - 0.49)
            #hiti = round(Int, hit)
            
            if hiti > 0
                phi = phi_interp_atan(y,x)
                theta = theta_interp[y,x]

                #rotatate star
                phi = phi - t*fs*2*pi
                phi = mod2pi(phi)
                
                img5[j,i] = painter(phi, theta)/2.0
                
                inside = spot(0.0, phi, theta)
            
                if inside
                    time = time_interp[y,x]

                    #compute 
                    earea = polyarea(x, y,
                                     frame_dxx, frame_dyy,
                                     phi, theta)
                
                    Xob = Xs_interp[y,x] 
                    cosa = cosa_interp[y,x]
                    dF, dE = radiation(Ir,
                                       x,y,
                                       phi, theta, cosa,
                                       X, Xob, Osb, sini, earea)
                    
                    #dF = flux_interp[y,x]
                    #dE = reds_interp[y,x]
                    
                    
                    #img4[j,i] = painter(phi, theta)
                    img5[j,i] += 3.0*dF/frame_dxdy
                    #img5[j,i] = 5.0
                    spot_flux[k] += dF
                end #inside spot
            end#hiti
        end #x
    end#y

    p10 = plot2d(img5, frame_xgrid, frame_ygrid, 0, 0, 5, "Blues")
    display(p10)

    
end#for t
toc()


#p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 0, "Blues")

