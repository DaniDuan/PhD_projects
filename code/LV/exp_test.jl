using LinearAlgebra, Random, Distributions, DifferentialEquations, Plots, CairoMakie, CSV, DataFrames, Turing, StatsPlots

include("./Model.jl")
# importing datasets
data_r = DataFrame(CSV.File("../../results/TPC/est_params_r.csv"))
data_a = DataFrame(CSV.File("../../results/TPC/est_params_a.csv"))

# setting reference temp to 10
r12 = data_r[:,:B0]
Ea_r = data_r[:,:Ea]
r0_all = temp_func((10+273.15), (12+273.15), r12, Ea_r)#, Ed_r, Th_r) 

row_Eα = data_a[:,:Ea]
Ea_α = reshape(repeat(row_Eα, 3),(3,3))
α0 = [-1.76997 2.94676 0.586733; -0.977813 -0.00855062 -0.649318;-0.621878 -1.25353 -0.777955]
# Plots.heatmap(α0, yflip=true, c = :gist_yarg,ticks = false)

# initial data
N = 3
Tr = 10+273.15
C0 = fill(0.01, N)
tspan = (0,24*24);
l = @layout [a b c]

for i in 1:3
    T = 10*i +273.15
    α = temp_func(T, Tr, α0, Ea_α)#, Ed_α, Th_α)
    r = temp_func(T, Tr, r0_all, Ea_r)#, Ed_r, Th_r)
    # α_12 = α[Not(1),Not(1)]
    # r_12 = r[Not(1)]
    # ## Simulation
    p = (N = N, α = α, r = r)
    prob = ODEProblem(GLV_model!, C0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))
    @eval $(Symbol("fig_temp$(10*i)")) = $(Plots.plot(sol, xaxis = "Time (hours)", yaxis = "Biomass", label=["S18" "W02" "W03"],
     legend=true, lc = [:darkgreen :blue :darkorange], title = string("Temp = ", i*10)))    
end
# color = darkgreen, blue, darkorange
Plots.plot(fig_temp10, fig_temp20, fig_temp30, layout = l, size = (1200, 400), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
savefig("../../results/Simulations/sim_All.png")


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
Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", label=["S18" "W02"], legend=true, lc = [:darkgreen :blue])    


Plots.scatter!(time, odedata, yerror=(ode_error, ode_error), label=["W02" "W03"], color = [:blue :darkorange],
 alpha = 0.7, markerstrokecolor=[:blue :darkorange])
