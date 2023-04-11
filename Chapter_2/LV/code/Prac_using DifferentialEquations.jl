using DifferentialEquations
using Plots


############# Example 1 - scalar equation
f(u,p,t) = 1.01*u
u0 = 0.5
tspan = (0.0, 1.0)

# solving
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# analyzing
sol[5]
sol.u
sol.t[8]
[t+u for (u,t) in tuples(sol)] # array of results by a comprehension over solution tuples
sol(0.45) # the solution at t = 0.45

# ploting
plot(sol, linewidth = 5, title = "solution to linear ODE",
 xaxis = "time", yaxis = "u(t)", label = "Line")
plot!(sol.t, t -> 0.5*exp(1.01t), lw = 3, ls=:dash, label = "true solution")


############### Example 2 - system

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
end

u0 = [1.0; 0 ; 0]
tspan = (0,100)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())

plot(sol, idxs = (1,2,3)) # plot on a 3D phase space
plot(sol, idxs = (0,2))


############### Parameterized functions !!!!!!!!!!!!

function parameterized_lorenz!(du, u, p, t)
    du[1] = p[1]*(u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end

u0 = [1.0, 0, 0]
tspan = (0,1)
p = [10.0, 28, 8/3]
prob = ODEProblem(parameterized_lorenz!, u0, tspan, p)
sol = solve(prob, Tsit5())


##################### Example 3


