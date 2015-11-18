using Winston


#image dimensions
mlim= 1.3
xmin = -1.8
xmax = mlim
ymin = -0.5
ymax = mlim


xlen = xmax-xmin
ylen = ymax-ymin
aspr=ylen/xlen

#variables
incl = pi/2.5
req = 1
u = 0.0
rg = u*req
muc = -rg/(req-rg)
  
R(theta) = 1.0
mu(phi,theta) = cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta)
b(phi,theta) = sqrt(req*(req+rg+req*mu(phi,theta)-rg*mu(phi,theta))/(1+mu(phi,theta)))

x(phi,theta) = b(phi,theta)*R(theta)*(cos(incl)*cos(theta)+cos(phi)*sin(incl)*sin(theta))
y(phi,theta) = b(phi,theta)*R(theta)*(sin(theta)*sin(phi))
z(phi,theta) = b(phi,theta)*R(theta)*(cos(theta)*sin(incl)-cos(incl)*cos(phi)*sin(theta))

p = FramedPlot()


#frontside
sopts = Dict()
sopts[:linekind] = "solid"
sopts[:color] = "black"

#backside
bsopts = Dict()
bsopts[:linekind] = "dotted"
bsopts[:color] = "black"

#second kind of lines (inclination & observer related)
isopts = Dict()
isopts[:linekind] = "dashed"
isopts[:color] = "black"

#axis style
axis_sopts = Dict()
axis_sopts[:linekind] = "solid"
axis_sopts[:color] = "black"
axis_sopts[:linewidth] = 3.5



function draw_longitude(p, phi, sopts;
                        start=0,
                        stop=pi/2,
                        rfac=1.0,
                        backside=false)
    xx = Float64[]
    yy = Float64[]
    for theta = linspace(start, stop, 200)
        if !backside
            if mu(phi,theta) >= muc
                push!(xx, y(phi,theta)*rfac)
                push!(yy, z(phi,theta)*rfac)
            else
                add(p, Curve(xx, yy, sopts))
                xx = Float64[]
                yy = Float64[]
            end
        else
            if mu(phi,theta) <= muc
                push!(xx, y(phi,theta)*rfac)
                push!(yy, z(phi,theta)*rfac)
            else
                add(p, Curve(xx, yy, sopts))
                xx = Float64[]
                yy = Float64[]
            end
        end
    end
    add(p, Curve(xx, yy, sopts))
    return p
end

function draw_latitude(p, theta, sopts; 
                       start=-pi,
                       stop=pi,
                       rfac=1.0,
                       backside=false)
    xx = Float64[]
    yy = Float64[]
    for phi = linspace(start, stop, 200)
        if !backside
            if mu(phi,theta) >= muc
                push!(xx, y(phi,theta)*rfac)
                push!(yy, z(phi,theta)*rfac)    
            else
                add(p, Curve(xx, yy, sopts))
                xx = Float64[]
                yy = Float64[]
            end
        else
            if mu(phi,theta) <= muc
                push!(xx, y(phi,theta)*rfac)
                push!(yy, z(phi,theta)*rfac)    
            else
                add(p, Curve(xx, yy, sopts))
                xx = Float64[]
                yy = Float64[]
            end
        end
    end
    add(p, Curve(xx, yy, sopts))
    return p
end

function draw_radial(p, phi, theta, sopts;
                     rfac=1.0)
    xx = Float64[0.0]
    yy = Float64[0.0]
    push!(xx, y(phi,theta)*rfac)
    push!(yy, z(phi,theta)*rfac)    
    add(p, Curve(xx, yy, sopts))

    return p
end

#draw circle around the sphere
#sopts2 = Dict()
#sopts2[:linekind] = "solid"
#sopts2[:color] = "red"

#t = linspace(0, 2pi, 500)
#xx = req*cos(t)
#yy = req*sin(t)
#add(p, Curve(xx, yy, sopts2))

#remove the axes
setattr(p, xrange=(xmin, xmax))
setattr(p, yrange=(ymin, ymax))
setattr(p.frame, 
        draw_spine=false, 
        draw_ticks=false,
        draw_ticklabels=false)


#coordinates of spot and observer
spot_theta=pi/5
spot_phi=0.67

obs_theta=incl
obs_phi=-pi/1.5

#helping curve for spot location
#p=draw_longitude(p, spot_phi-0.07, sopts)
#p=draw_longitude(p, spot_phi+0.03 , sopts)
p=draw_longitude(p, spot_phi , sopts)

#borders
p=draw_longitude(p, pi/2, sopts)
p=draw_longitude(p, -pi/2, sopts)

p=draw_latitude(p, pi/2, sopts)
p=draw_latitude(p, pi/2, bsopts, backside=true)

#spot position
p=draw_radial(p, spot_phi, spot_theta, sopts)
p=draw_longitude(p, spot_phi, sopts, start=0, stop=spot_theta, rfac=0.17) #theta angle
p=draw_longitude(p, spot_phi, sopts, start=0, stop=spot_theta, rfac=0.17, backside=true) #theta angle
#add(p, PlotLabel(.53, .58, "<i>\\theta</i>"))


#xyz axis
p=draw_radial(p, obs_phi, pi/2, axis_sopts, rfac=1.4) #x
p=draw_radial(p, obs_phi+pi/2, pi/2, axis_sopts, rfac=1.4) #y
p=draw_radial(p, 0.0, 0.0, axis_sopts, rfac=1.4) #z

#help axis
p=draw_radial(p, spot_phi, pi/2, sopts)
#p=draw_radial(p, 0.0, pi/2, sopts)


#observer
#p=draw_radial(p, obs_phi, obs_theta, isopts, rfac=1.0)
p=draw_radial(p, obs_phi, obs_theta, isopts, rfac=1.5)
p=draw_radial(p, obs_phi, pi/2, isopts)
#p=draw_latitude(p, obs_theta, isopts)
#p=draw_latitude(p, obs_theta, isopts, backside=true)
p=draw_longitude(p, obs_phi, sopts, start=0, stop=obs_theta, rfac=0.1) #incl angle
p=draw_longitude(p, obs_phi, sopts, start=0, stop=obs_theta, rfac=0.1, backside=true) #incl angle
#add(p, PlotLabel(.47, .55, "<i>i</i>"))

p=draw_latitude(p, pi/2, sopts, start=spot_phi, stop=obs_phi, rfac=0.15) #phi angle
p=draw_latitude(p, pi/2, sopts, start=spot_phi, stop=obs_phi, rfac=0.15, backside=true) #phi angle
#add(p, PlotLabel(.47, .45, "<i>\\phi</i>"))

setattr(p, aspect_ratio=aspr)
Winston.file(p, "fig_geometry.eps")
