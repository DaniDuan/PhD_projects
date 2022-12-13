using DifferentialEquations
using Plots
# using RCall # for fun

include("./Model.jl")
# test
N = 2
T = 0+273.15
p = Params(N, T)
C0 = fill(0.1, N)
tspan = (1,1000)


# solving test
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())
sol.u

Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false)
# savefig("./test.png")

# # for fun
# $
# library(beepr)
# beep(sound = 4, expr = NULL) 


################## 
N = 2
C0 = fill(0.1, N)
tspan = (1,1000)
Temps = [0, 10, 25]
fig = Figure(resolution = (1200, 400))
for i in 1:3
    T = Temps[i]+273.15
    Random.seed!(1)    
    p = Params(N, T)
    prob = ODEProblem(GLV_model!, C0, tspan, p);
    sol = solve(prob, Tsit5());
    max_value = maximum(maximum(sol.u))
    odeSol(x, y) = Point2f(p.r[1]+p.α[1,1]*x+ p.α[1,2]*y, p.r[2]+p.α[2,2]*y+ p.α[2,1]*x) # Point2f(x',y')
    ax = Axis(fig[1, i], xlabel = "C1", ylabel = "C2", backgroundcolor = :white)
    streamplot!(ax, odeSol, 0 .. max_value, 0 .. max_value, colormap = Reverse(:Blues), gridsize = (32, 32), arrow_size = 10)
end

# Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false)    

save("./test.png", fig)
