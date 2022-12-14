using CairoMakie
include("./temp_params.jl")

N = 2
T = 0+273.15
p = Params(N, T)

function phase_plane_for2(p, max_value)
    odeSol(x, y) = Point2f(p.r[1]+p.α[1,1]*x+ p.α[1,2]*y, p.r[2]+p.α[2,2]*y+ p.α[2,1]*x) # Point2f(x',y')
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], xlabel = "C1", ylabel = "C2", backgroundcolor = :white)
    streamplot!(ax, odeSol, 0 .. max_value*1.5, 0 .. max_value*1.5, colormap = Reverse(:Blues),
        gridsize = (32, 32), arrow_size = 10)
    return (fig)
end

# fig = phase_plane_for2(p, 2.5)
# save("phase_test_0.png", fig)
