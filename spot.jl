#Compute spot on the NS


#Load photon bender
include("bender.jl")

#Compute raw image
#include("comp_img.jl")

#Interpolate from raw image
#include("img.jl")



########
function spot(t, phi, theta;
              stheta = pi/3, #spot latitude
              delta = 40*pi/180 #spot half-opening angle
              )

    #Vincenty's formula
    d = great_circle_dist(0.0, phi, stheta, theta)
    
    if abs(d) < delta
        return true
    end
    
    return false
end

img4 = zeros(Ny_dense, Nx_dense) #debug array
Ir(cosa) = 1.0 #isotropic beaming
Nt = 60
times = linspace(0, 1/fs, Nt)


tic()
for k = 1:Nt
    img4[:,:] = 0.0
    
    t = times[k]
    for j = y1s:y2s
        y = y_grid_d[j]
        for i = x1s[j]:x2s[j]
            x = x_grid_d[i]

            theta = theta_interp[y,x]

            #rotate star
            phi = phi_interp_atan(y,x)
            phi = phi - t*fs*2*pi
            phi = mod2pi(phi)
      

            img4[j,i] = painter(phi, theta)/10

            inside = spot(0.0, phi, theta)
            
            if inside
            #if true    
                time = time_interp[y,x]
                Xob = Xs_interp[y,x] 
                cosa = cosa_interp[y,x]
                dF, dE = radiation(Ir,
                                   x,y,
                                   phi, theta, cosa,
                                   X, Xob, Osb, sini, img3[j,i])

                
                #img4[j,i] = painter(phi, theta)
                img4[j,i] += dF*1e4
                
            end #if inside            
        end# for x
    end#for y

    p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 5, "Blues")
    #p10 = plot2d(img4, x_grid_d, y_grid_d, 0, -1., 1., "RdBu")
    display(p10)
    
end#for t
toc()



p10 = plot2d(img4, x_grid_d, y_grid_d, 0, 0, 0, "Blues")

