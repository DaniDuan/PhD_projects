using LinearAlgebra, Random, Distributions, DifferentialEquations, Plots, CairoMakie, CSV, DataFrames, Turing, StatsPlots

include("./Model.jl")
# importing datasets
data_r = DataFrame(CSV.File("../../results/TPC/est_params_r.csv"))
data_a = DataFrame(CSV.File("../../results/TPC/est_params_a.csv"))
data_ref = DataFrame(CSV.File("../../results/TPC/ref_biomass.csv"))

# setting reference temp to 10
r12 = data_r[:,:B0]
Ea_r = data_r[:,:Ea]
Ed_r = fill(4, 3)
Th_r = data_r[1:3,:Th] .+ 273.15 # from data
r0_all = temp_func((10+273.15), r12, Ea_r, Ed_r, Th_r) 

α0_v = -data_a[:,:B0]
row_Eα = data_a[:,:Ea]
Ed_α = fill(4, 3, 3)
row_Th = data_a[:,:Th]
α0_all = temp_func((10+273.15), α0_v, row_Eα, Ed_α, row_Th)

# data_r = data_r[Not(2), :]
# data_a = data_a[Not(2), :]

######## S18 + W03 ########
N = 2
C0 = fill(0.01, N)
tspan = (0,100)
T = 10+273.15

# r0 = r0_all[1:N]
r0 = r0_all[Not(1)]
# α0ii = α0_all[1:N]
α0ii = α0_all[Not(1)]

# #----------interspecies interaction is random for now---------------#
# α0_range = -α0ii*0.25
# α0 = reduce(vcat, transpose.([rand(Normal(α0ii[i], α0_range[i]), N) for i in 1:N]))
# #----------interspecies interaction is random for now---------------#
# α0[diagind(α0)] = α0ii
# Ea_α = reshape(repeat(row_Eα, N),(N,N))
# Th_α = [row_Th row_Th row_Th] .+ 273.15
# α = temp_func(T, α0, Ea_α, Ed_α, Th_α)

# ## Simulation
# p = (N = N, α = α, r = r)
# prob = ODEProblem(GLV_model!, C0, tspan, p)
# sol = solve(prob, Tsit5())

# ## Analysis
# δT = -1/0.0000862 * (1/T - 1/(12+273.15))
# Δr0 = log(r0[1]/r0[2])
# Δα0 = log.(α0[1,:]./α0[2,:])
# Ei = Ea_r[1] - row_Eα[1]
# Ej = Ea_r[2] - row_Eα[2]
# fst = (Ei - Ej)*δT - (Δα0[2] - Δr0)
# snd = Δα0[1] - Δr0 - (Ei - Ej)*δT

# print("richness:", count(x-> x > 1e-5, sol.u[length(sol)]), 
# "\nbiomass: ", sol.u[length(sol)], "\n", fst > 0, " ", snd > 0)

# color = ("darkgreen", "blue", "chocolate2")
# Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", label=["S18" "W03"], legend=true, lc = [:darkgreen :darkorange])    

function LV(du, u, p, t)
    # Model parameters.
    α_ij, α_ji = p
    # Current state.
    x, y = u

    # Evaluate differential equations.
    du[1] = x*(r0[1]+α0ii[1]*x + α_ij*y)
    du[2] = y*(r0[2]+α0ii[2]*x + α_ji*y)
end

odedata = hcat(data_ref[:,2], data_ref[:,3])

time = [1, 24, 49, 72, 91]

α_ij = -0.5
α_ji = -1

α_ij/α0ii[2] - r0[1]/r0[2] #<0
α0ii[1]/α_ji - r0[1]/r0[2] #>0

p = [α_ij α_ji]
prob = ODEProblem(LV, C0, tspan, p)
sol = solve(prob, Tsit5(); saveat=0.5)
Plots.plot(sol)
Plots.scatter!(time, odedata)

@model function fitlv(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α_ij ~ truncated(Normal(-0.5, 0.5); lower=-1, upper=0)
    α_ji ~ truncated(Normal(-1, 0.5); lower=-1.5, upper=-0.5)

    # Simulate Lotka-Volterra model. 
    p = [α_ij,α_ji]
    predicted = solve(prob, Tsit5(); p=p, saveat=0.5)

    # Observations.
    if length(predicted) == 201
        for i in 1:5
            t = time[i]
            data[i,:] ~ MvNormal(predicted[t], σ^2 * I)
        end
    end
    return nothing
end

# data = Matrix{Union{Missing, Float64}}(missing, 101, 2)
# for i in 1:5
#     t = time[i]
#     data[t,:] = odedata[i,:]
# end

model = fitlv(odedata, prob)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

Plots.plot(chain)
