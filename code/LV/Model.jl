using DifferentialEquations
using Plots

include("./temp_params.jl")

# test parameters
N = 2;
r = [0.1, 0.4];
α = [-0.2 -0.1;
     -0.1 -0.2]
tspan = (0,100);
C0 = fill(0.1, N)

p = (N = N, r = r, α = α)

# the ODE function
function GLV_model!(dx, x, p, t)
    dx .= x.*(p.r+p.α*x)
end

# solving test
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())
