using Winston
using Colors
using Grid
using Interpolations

include("plot2d.jl")

radius = 3.0
grid_x_size = 7.0
grid_y_size = 7.0
grid_nx = 100
grid_ny = 100
grid_dx = grid_x_size/grid_nx
grid_dy = grid_y_size/grid_ny
xmax = grid_x_size/2
xmin = -xmax
ymax = grid_y_size/2
ymin = -ymax

x_grid = linspace(xmin, xmax, grid_nx)
y_grid = linspace(ymin, ymax, grid_ny)'
xy_grid = [(x,y) for x in x_grid, y in y_grid] 
xx_grid = repmat(x_grid', grid_nx, 1)
yy_grid = repmat(y_grid'[end:-1:1], 1, grid_ny)

#println("xy_grid $xy_grid")

# check radial distance of points from center against sphere size
rdist_grid = sqrt(x_grid.^2 .+ y_grid.^2)
inside  = find(rdist_grid .<= radius)
outside = find(rdist_grid .> radius)
#println("size(rdist_grid)", size(rdist_grid))

# get corresponding phis and thetas
zz_grid = zeros(rdist_grid)
zz_grid[inside] = sqrt(radius^2 .- rdist_grid[inside].^2)

theta_grid = acos(zz_grid ./ sqrt(xx_grid.^2 .+ yy_grid.^2 .+ zz_grid.^2))
phi_grid = atan2(yy_grid, xx_grid)
theta_grid[outside] = pi/2
phi_grid[outside] = 0

plot_z     = plot2d(zz_grid    , x_grid , y_grid)
plot_theta = plot2d(theta_grid , x_grid , y_grid)
plot_phi   = plot2d(phi_grid   , x_grid , y_grid)

# try lambert azimuthal equal area projection to get the differential areas
function lambert(phis, thetas)
    Rs = 2 * cos((pi-thetas) ./ 2)
    xs = cos(phis) .* Rs
    ys = sin(phis) .* Rs
    xys = collect(zip(xs, ys))
    xys2 = append!(xys[2:end], [xys[1]])
    polygon = collect(zip(xys, xys2))
    area = 0.5 * abs(sum(
        [x0*y1 - x1*y0 for ((x0,y0), (x1,y1)) in polygon]))
    return area
end

function lambert2(xs, ys, zs)
    xs /= radius
    ys /= radius
    zs /= radius
    zs = -zs
    Xs = sqrt(2 ./(1 .- zs)) .* xs
    Ys = sqrt(2 ./(1 .- zs)) .* ys
    XYs = collect(zip(Xs, Ys))
    XYs2 = append!(XYs[2:end], [XYs[1]])
    polygon = collect(zip(XYs, XYs2))
    area = 0.5 * abs(sum(
        [x0*y1 - x1*y0 for ((x0,y0), (x1,y1)) in polygon]))
    return area
end

function lambert3(phis, thetas, phi0, theta0)
    lphis = pi/2 - thetas
    lphi0 = pi/2 - theta0
    llambdas = phis
    llambda0 = phi0

    kprime = sqrt(2 ./ (1 .+ sin(lphi0) .* sin(lphis) .+ cos(lphi0) .* cos(lphis) .* cos(llambdas-llambda0)))
    xs = kprime .* cos(lphis) .* sin(llambdas-llambda0)
    ys = kprime .* (cos(lphi0) .* sin(lphis) .- sin(lphi0) .* cos(lphis) .* cos(llambdas-llambda0))

    xys = collect(zip(xs, ys))
    xys2 = append!(xys[2:end], [xys[1]])
    polygon = collect(zip(xys, xys2))
    area = 0.5 * abs(sum(
        [x0*y1 - x1*y0 for ((x0,y0), (x1,y1)) in polygon]))
    return area
end

# compute areas on sphere for each point inside the circle
area_grid = zeros(xx_grid)
area_grid2 = zeros(xx_grid)
area_grid3 = zeros(xx_grid)
for index in inside
    (xi, yi) = ind2sub(xx_grid, index)
    # get a quad from the grid
    p1 = (xi-1, yi-1)
    p2 = (xi+1, yi-1)
    p3 = (xi+1, yi+1)
    p4 = (xi-1, yi+1)
    pts = [p1, p2, p3, p4]

    # check that all points are inside
    pts = filter(x -> sub2ind(size(xx_grid), x[1], x[2]) in inside, pts)

    # convert to phis and thetas
    phis = Float64[phi_grid[xind,yind] for (xind,yind) in pts]
    thetas = Float64[theta_grid[xind,yind] for (xind,yind) in pts]
    xs = Float64[xx_grid[xind,yind] for (xind,yind) in pts]
    ys = Float64[yy_grid[xind,yind] for (xind,yind) in pts]
    zs = Float64[zz_grid[xind,yind] for (xind,yind) in pts]
    area_grid[index] = lambert(phis, thetas)
    area_grid2[index] = lambert2(xs, ys, zs)
    area_grid3[index] = lambert3(phis, thetas, 0, 0)
end
plot_area  = plot2d(area_grid,  x_grid , y_grid)
plot_area2 = plot2d(area_grid2, x_grid , y_grid)
plot_area3 = plot2d(area_grid3, x_grid , y_grid)


# plot everything
#plot_array = [plot_z, plot_theta, plot_phi, plot_area, plot_area2]
plot_array = [plot_area, plot_area2, plot_area3]
ps = size(plot_array, 1)
println("plot array size $ps")
plot_table = Table(1,ps)
for pair in enumerate(plot_array)
    println("pair $pair pair[1] ", pair[1], " pair[2] ", pair[2])
    plot_table[1,pair[1]] = pair[2]
end
display(plot_table)


