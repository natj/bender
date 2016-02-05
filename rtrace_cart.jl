#Compute image using 4 quadrants;
#(can be parallellized easily)

######################
#grid setup
xmin = -10
xmax = -xmin +0.005
ymin = -10
ymax = -ymin +0.005

Ny = 100
Nx = 100

xmid = round(Int,Nx/2)
ymid = round(Int,Ny/2)
hdata3 = zeros(Ny, Nx)

#dx = (xmax-xmin)/(Nx-1)
#dy = (ymax-ymin)/(Ny-1)
#x_grid = collect(xmin:dx:xmax)
#y_grid = collect(ymin:dy:ymax)   

x_grid = collect(linspace(xmin,xmax, Nx))
y_grid = collect(linspace(ymin,ymax, Ny))
dx = diff(x_grid)[1]
dy = diff(y_grid)[1]


Times = zeros(Ny,Nx)
Phis = zeros(Ny,Nx)
Thetas = zeros(Ny,Nx)
hits = zeros(Ny,Nx)
cosas = zeros(Ny,Nx)
Xs = zeros(Ny,Nx)

#Get middle point
time, phi, theta, Xob, hit, cosa = bender3(x_grid[xmid], y_grid[ymid], sini,
                                           X, Osb,
                                           beta, quad, wp, Rg)

Times[ymid, xmid] = time
Phis[ymid, xmid] = phi
Thetas[ymid, xmid] = theta


println("Computing image...")
tic()
#Upper region
for j = ymid:Ny
    #upper right corner
    if mod(j, 20) == 0; println("u: $j"); end

    for i = xmid:Nx
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa
        if !hit
            break
        end
    end

    #upper left corner
    for i = xmid:-1:1
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        if !hit
            break
        end
    end
end

#Lower region
for j = ymid:-1:1
    if mod(j, 20) == 0; println("l: $j"); end
    #lower right corner
    for i = xmid:Nx
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                             X, Osb,
                                             beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        if !hit
            break
        end
    end

    #lower left corner
    for i = xmid:-1:1
        time, phi, theta, Xob, hit, cosa = bender3(x_grid[i], y_grid[j], sini,
                                                   X, Osb,
                                                   beta, quad, wp, Rg)
        Times[j,i] = time
        Phis[j,i] = phi
        Thetas[j,i] = theta
        Xs[j,i] = Xob
        hits[j,i] = float(hit)
        cosas[j,i] = cosa

        if !hit
            break
        end
    end
end
toc()
