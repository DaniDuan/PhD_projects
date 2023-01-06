using LinearAlgebra, Random, Distributions, DifferentialEquations, Plots, CairoMakie, CSV, DataFrames
include("./Model.jl")
data_r = DataFrame(CSV.File("../../results/TPC/est_params_r.csv"))
data_a = DataFrame(CSV.File("../../results/TPC/est_params_a.csv"))

data_r = data_r[Not(2), :]
data_a = data_a[Not(2), :]

######## S18 + W02 ########
N = 2
C0 = fill(0.01, N)
tspan = (1,25)
T = 20+273.15

r0 = data_r[1:N,:B0]
Ea_r = data_r[1:N,:Ea]
Ed_r = fill(4, N)
Th_r = data_r[1:N,:Th] .+ 273.15 # from data
r = temp_func(T, r0, Ea_r, Ed_r, Th_r) 
α0_v = -data_a[1:N,:B0]
#----------interspecies interaction is random for now---------------#
α0_range = data_a[1:N,:B0]*0.25
α0 = reduce(vcat, transpose.([rand(Normal(α0_v[i], α0_range[i]), N) for i in 1:N]))
#----------interspecies interaction is random for now---------------#
α0[diagind(α0)] = -data_a[1:N,:B0]
row_Eα = data_a[1:N,:Ea]
Ea_α = reshape(repeat(row_Eα, N),(N,N))
Ed_α = fill(4, N, N)
row_Th = data_a[1:N,:Th]
Th_α = [row_Th row_Th row_Th] .+ 273.15
α = temp_func(T, α0, Ea_α, Ed_α, Th_α)

## Simulation
p = (N = N, α = α, r = r)
prob = ODEProblem(GLV_model!, C0, tspan, p)
sol = solve(prob, Tsit5())

## Analysis
δT = -1/0.0000862 * (1/T - 1/(12+273.15))
Δr0 = log(r0[1]/r0[2])
Δα0 = log.(α0[1,:]./α0[2,:])
Ei = Ea_r[1] - row_Eα[1]
Ej = Ea_r[2] - row_Eα[2]
fst = (Ei - Ej)*δT - (Δα0[2] - Δr0)
snd = Δα0[1] - Δr0 - (Ei - Ej)*δT

print("richness:", count(x-> x > 1e-5, sol.u[length(sol)]), 
"\nbiomass: ", sol.u[length(sol)], "\n", fst > 0, " ", snd > 0)

color = ("darkgreen", "blue", "chocolate2")
Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", label=["S18" "W03"], legend=true, lc = [:darkgreen :darkorange])    
