using LinearAlgebra, Random, Distributions, DifferentialEquations, Plots, CairoMakie
include("./Model.jl")

tspan = (1,2500)
C0 = fill(0.1, 2)
det_value = -[-0.01 0.00 0.01]
l = @layout [a b c]
fig = Figure(resolution = (1200, 400))

for i in 1:3
    r = rand(Normal(0.21, 0.01), 2) # from data
    α = rand(Normal(-1.18, 0.5), 2,2)
    α[2,1] = (det(diagm(diag(α)))+det_value[i])/α[1,2]
    if i != 2
        while α[1,2]/α[2,2] >= r[1]/r[2] || α[1,1]/α[2,1] <= r[1]/r[2] # α[1,2]/α[2,2] < r[1]/r[2] and α[1,1]/α[2,1] > r[1]/r[2]
        r = rand(Normal(0.21, 0.01), 2) # from data
        α = rand(Normal(-1.18, 0.5), 2,2)
        α[2,1] = (det(diagm(diag(α)))+det_value[i])/α[1,2]
        end
    end
    print(det(α), "\n")
    p = (N = 2, α = α, r = r)
    prob = ODEProblem(GLV_model!, C0, tspan, p);
    sol = solve(prob, Tsit5());
    @eval $(Symbol("fig_dyn$i")) = $(Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false))
    max_value = maximum(maximum(sol.u))
    odeSol(x, y) = Point2f(p.r[1]+p.α[1,1]*x+ p.α[1,2]*y, p.r[2]+p.α[2,2]*y+ p.α[2,1]*x) # Point2f(x',y')
    ax = Axis(fig[1, i], xlabel = "C1", ylabel = "C2", backgroundcolor = :white)
    streamplot!(ax, odeSol, 0 .. max_value*1.5, 0 .. max_value*1.5, colormap = Reverse(:Blues), gridsize = (32, 32), arrow_size = 10)
end
display(Plots.plot(fig_dyn1, fig_dyn2, fig_dyn3, layout = l, size = (950, 300)))
display(fig)