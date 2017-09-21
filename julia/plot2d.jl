function plot2d(sdata, xx_grid, yy_grid,
                smooth=0,
                hmin=0.0,
                hmax=0.0,
                colmap="RdBu";
                xlabel="",
                ylabel="")
    
    xmin = minimum(xx_grid)
    xmax = maximum(xx_grid)
    ymin = minimum(yy_grid)
    ymax = maximum(yy_grid)
    
    p = FramedPlot(xlabel=xlabel,
                   ylabel=ylabel)
    #hmin = minimum(sdata[find(sdata)])

    if hmin == 0 && hmax == 0
        hmin = minimum(sdata)
        hmax = maximum(sdata)
    end
    
    println("min=$hmin max=$hmax")
    clims = (hmin, hmax)

    #colormap
    #cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Colors.colormap(colmap))]
    #cm = Uint32[Color.convert(Color.RGB24,c) for c in flipud(Colors.colormap("Blues"))]
    #cm = UInt32[Colors.convert(Colors.RGB24,c) for c in Colors.colormap(colmap)]
    cm = UInt32[Colors.convert(Colors.RGB24,c) for c in flipdim(Colors.colormap(colmap),1)]
    unshift!(cm, 0x00ffffff)
    #push!(cm, 0x00ffffff)

    # 5x5 gaussian
    kernel2d=(1.0/273.)*[1.0 4.0 7.0 4.0 1.0;
                             4.0 16. 26. 16. 4.0;
                             7.0 26. 41. 26. 7.0;
                             1.0 4.0 7.0 4.0 1.0;
                             4.0 16. 26. 16. 4.0]                 

    kernelx, kernely = size(kernel2d)
    kx = kernelx-3
    ky = kernely-3
    
    for i in 1:smooth
        sdata = conv2(sdata*1.0, kernel2d)
        sdata = sdata[ky+1:end-ky, kx+1:end-ky] #resize data to remove edges
    end
    
    #image
    img = Winston.data2rgb(sdata, clims, cm)
    add(p, Image((xmin, xmax), (ymin, ymax), img;))
    setattr(p, xrange=(xmin, xmax))
    setattr(p, yrange=(ymin, ymax))
    #setattr(p, xlabel=xlabel)
    #setattr(p, ylabel=ylabel)
        
    return p
end
#########################


function chess_board(phi, theta, cosa=0, dF=0, Es=0)

    none = 0
    white = 1
    black = 2
    
    #if abs(theta) < pi/20
    #    return black
    #end
    #if !(0 < theta < pi)
    #    return black
    #end

    #println("phi = $phi")
    #println("theta = $theta")
    #xd = round(Int,80*phi/2pi)
    #yd = round(Int,80*theta/2pi)
    xd = round(Int,60*phi/2pi)
    yd = round(Int,60*theta/2pi)

    if isodd(xd) && isodd(yd)
        return black
    elseif iseven(xd) && iseven(yd)
        return black
    else
        return white
    end
end

function polar_caps(phi, theta, cosa=0, dF=0, Es=0)
    none = 0
    white = 1
    black = 2
    
    size = pi/10
    if (-size < theta < size)
        return black
    elseif (pi-size < theta < pi+size) 
        return black
    else
        return white
    end
end

