#test fucked up polar angles

N = 9

r = 2.0
chis = linspace(0.0, 2pi, 9)

xs = zeros(N)
ys = zeros(N)

for i = 1:N
    x = r*sin(chis[i])
    y = r*cos(chis[i])

    println("x:$x y:$y")
    xs[i] = x
    ys[i] = y
end

println()

for i=1:N

    x = xs[i]
    y = ys[i]

    r = hypot(x,y)
    chi = mod2pi(pi/2 - atan2(y,x))

    println("x:$x y:$y r=$r chi=$(chi/pi)")
end


phi0 = -2.9123889170923705
the0 = 0.42946959557816783
sinphis = [0.0,0.42246406473560055,0.3649600332965458,0.0]
cosphis = [1.0,1.0,1.0,1.0]
thetas = [-0.6040411556929683,-2.9302058807662563,-3.431687364793615,-0.7089982632590643]
area_sphere_lambert2(phi0, the0, sinphis, cosphis, thetas)

println("sxxxx")

phi0 = -2.8913935542943707
the0 = 0.38023520468294647
sinphis = [0.0,0.41805020983662955,0.3610066075301765,0.0]
cosphis = [1.0,1.0,1.0,1.0]
thetas = [-0.6063085387422033,-2.9238509673642397,-3.4259478665419705,-0.7120003938552091]

area_sphere_lambert2(phi0, the0, sinphis, cosphis, thetas)
#[1.016276183492174,11.152865154640086,-1359.0032960795247,1.032035037298757]
