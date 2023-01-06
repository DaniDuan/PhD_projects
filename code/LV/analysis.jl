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

############################################################################################
N = 2
C0 = fill(0.1, N)
tspan = (1,5000)
Ed_r = fill(4, N)
Th_r = rand(Normal(25+273.15, 5), N) # from data
Ed_α = fill(4, N, N)
Th_α = rand(Normal(24+273.15, 1), N, N) # test variation

# all B0s 
r0 = rand(Normal(0.21, 0.01), N) # from data
α0 =  0
while det(α0) <= 0
    α0 =  rand(Normal(-1.18, 0.5), N,N) # mean value from data, test variation
end

Δr0 = log(r0[1]/r0[2])
Δα0 = log.(α0[1,:]./α0[2,:])


T = 20+273.15
δT = -1/0.0000862 * (1/T - 1/(12+273.15))

# all Eas
Ea_r = rand(Normal(0.95, 0.16), N) # from data
ran_Eα = rand(Normal(2.2, 0.5), N)
Ea_α = [ran_Eα ran_Eα] # mean value from data, test variation

Ei = Ea_r[1] - ran_Eα[1]
Ej = Ea_r[2] - ran_Eα[2]

fst = (Ei - Ej)*δT - (Δα0[2] - Δr0)
snd = Δα0[1] - Δr0 - (Ei - Ej)*δT

# running the model
r = temp_func(T, r0, Ea_r, Ed_r, Th_r) 
α = temp_func(T, α0, Ea_α, Ed_α, Th_α)
p = (N = N, α = α, r = r)
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())

print("richness:", count(x-> x > 1e-5, sol.u[length(sol)]), 
"\nbiomass: ", sol.u[length(sol)], "\n", fst > 0, " ", snd > 0)
