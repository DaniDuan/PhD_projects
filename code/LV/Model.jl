using DifferentialEquations
using Plots
using Random, Distributions

include("./temp_params.jl")

# test
N = 5
T = 15+273.15
p = Params(N, T)

# the ODE function
function GLV_model!(dx, x, p, t)
    dx .= x.*(p.r+p.α*x)
end

# solving test
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())
