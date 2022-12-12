using DifferentialEquations
using Plots
using RCall # for fun

include("./temp_params.jl")
# test
N = 2
T = 15+273.15
p = Params(N, T)
C0 = fill(0.5, N)
tspan = (1,100)


# solving test
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())
sol.u

plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false)
# savefig("./test.png")

# for fun
$
library(beepr)
beep(sound = 4, expr = NULL) 