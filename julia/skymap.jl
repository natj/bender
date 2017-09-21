
incls  = Float64[5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
colats = Float64[10.0, 30.0, 50.0, 70.0, 90.0]

for incl_deg in reverse(incls)

    println("incl: ",incl_deg)
    
    incl = deg2rad(incl_deg)
    include("bender.jl")
    include("rtrace.jl")
    include("radiation.jl")
    include("img.jl")

    for colat_deg in colats

        println("colat: ",colat_deg)
        
        colat = deg2rad(colat_deg)
        include("spot.jl")
        include("spot2nd.jl")
    end
end
