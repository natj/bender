#Compute and save trace of photons
#Used for 3d_photon_path.py to plot the paths

#grid setup
#x_grid = collect(linspace(-10, 10, 10))
#y_grid = collect(linspace(-10, 10, 10))

x_grid = collect(linspace(-10.0, 10, 200))
y_grid = [4.0]
    
mkpath("ppath")
    
for j = 1:length(y_grid)
    ypoint = y_grid[j]
    #ypoint = 0.005
    
    for i = 1:length(x_grid)
        xpoint = x_grid[i]
        #xpoint = 4.998743743743744
            
        rns, yns, zns, ers, lvs, hit = bender3(xpoint, ypoint, sini,X, Osb,
                                               beta, quad, wp, Rg)

        if hit == 1.0
            fname = "ppath/p_$(round(xpoint,3))_$(round(ypoint,3)).csv"
            println(fname)
            wmatr = hcat(rns, yns, zns, ers, lvs)
            writecsv(fname, wmatr)
        end
    end
end
