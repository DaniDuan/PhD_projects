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
    sol = solve(prob, AutoTsit5(Rosenbrock23()));
    @eval $(Symbol("fig_dyn$i")) = $(Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false))
    max_value = maximum(maximum(sol.u))
    odeSol(x, y) = Point2f(p.r[1]+p.α[1,1]*x+ p.α[1,2]*y, p.r[2]+p.α[2,2]*y+ p.α[2,1]*x) # Point2f(x',y')
    ax = Axis(fig[1, i], xlabel = "C1", ylabel = "C2", backgroundcolor = :white)
    streamplot!(ax, odeSol, 0 .. max_value*1.5, 0 .. max_value*1.5, colormap = Reverse(:Blues), gridsize = (32, 32), arrow_size = 10)
end
display(Plots.plot(fig_dyn1, fig_dyn2, fig_dyn3, layout = l, size = (950, 300)))
display(fig)

############################################################################################
## coexistence at all temps
N = 2
C0 = fill(0.01, N)
tspan = (0,1000)
Temps = [0, 12, 25]
Tr = 12+273.15
l = @layout [a b c]
fig = Figure(resolution = (1200, 400))
# δT = -1/0.0000862 * (1/T - 1/(10+273.15))
fst, snd, r0, α0, Ea_r, Ea_α, Ei, Ej = zeros(8)
while fst<= 0 || snd<=0 
    # all B0s 
    r0 = rand(Normal(0.17, 0.005), N) # from data
    α0 =  rand(Normal(-0.28, 1.91), N,N) # from data
    α0_diag = rand(truncated(Normal(-0.85, 0.78); upper = 0), N)
    α0[diagind(α0)] = α0_diag
    while det(α0)<=0
        α0 =  rand(Normal(-0.28, 1.91), N,N) # from data
        α0_diag = rand(truncated(Normal(-0.85, 0.78); upper = 0), N)
        α0[diagind(α0)] = α0_diag    
    end
    # all Eas
    Ea_r = rand(Normal(0.95, 0.17), N) # from data
    ran_Eα = rand(truncated(Normal(2.19, 3.74); lower = 0, upper = 5), N) # from data
    Ea_α = reshape(repeat(ran_Eα,N), N,N)
    Ei = Ea_r[1] - ran_Eα[1]
    Ej = Ea_r[2] - ran_Eα[2]
    # fst = exp((Ei - Ej)*δT)- r0[2]*α0[1,2]/(r0[1]*α0[2,2])
    # snd = exp((Ej - Ei)*δT)- r0[1]*α0[2,1]/(r0[2]*α0[1,1])
    if isless(Ei,Ej) T_ex1,T_ex2 = 25,0 else T_ex1,T_ex2 = 0,25 end
    fst = exp((Ei - Ej)*(-1/0.0000862 * (1/(T_ex1+273.15) - 1/(10+273.15)))) - r0[2]*α0[1,2]/(r0[1]*α0[2,2])
    snd = exp((Ej - Ei)*(-1/0.0000862 * (1/(T_ex2+273.15) - 1/(10+273.15)))) - r0[1]*α0[2,1]/(r0[2]*α0[1,1])
end

for i in 1:3
    T = Temps[i]+273.15
    δT = -1/0.0000862 * (1/T - 1/(10+273.15))
    # running the model
    r = temp_func(T, Tr, r0, Ea_r)
    α = temp_func(T, Tr, α0, Ea_α)
    p = (N = N, α = α, r = r)
    prob = ODEProblem(GLV_model!, C0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))
    fst = exp((Ei - Ej)*δT) - r0[2]*α0[1,2]/(r0[1]*α0[2,2])
    snd = exp((Ej - Ei)*δT) - r0[1]*α0[2,1]/(r0[2]*α0[1,1])
    print("richness:", count(x-> x > eps(), sol.u[length(sol)]), "\nbiomass: ", sol.u[length(sol)], "\n", fst > 0, " ", snd > 0, "\n", α, "\n")
    @eval $(Symbol("fig_dyn$i")) = $(Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", label=["i" "j"], legend=true, title = string("Temp = ", Temps[i])))
    max_value = maximum(maximum(sol.u))
    odeSol(x, y) = Point2f(p.r[1]+p.α[1,1]*x+ p.α[1,2]*y, p.r[2]+p.α[2,2]*y+ p.α[2,1]*x) # Point2f(x',y')
    ax = CairoMakie.Axis(fig[1, i], xlabel = "C1", ylabel = "C2", backgroundcolor = :white)
    streamplot!(ax, odeSol, 0 .. max_value*1.5, 0 .. max_value*1.5, colormap = Reverse(:Blues), gridsize = (32, 32), arrow_size = 10)
end
# CairoMakie.save("../../results/Simulations/coexistence_phase.png", fig)
Plots.plot(fig_dyn1, fig_dyn2, fig_dyn3, layout = l, size = (1200, 400), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
# savefig("../../results/Simulations/coexistence_biomass.png")
# display(fig)

print("richness:", count(x-> x > 1e-5, sol.u[length(sol)]), 
"\nbiomass: ", sol.u[length(sol)], "\n", fst > 0, " ", snd > 0, "\n", α)

# Plots.plot(sol)

######################### !!! Coexistence Condition #####################################

using LinearAlgebra, Random, Distributions, DifferentialEquations, Plots, CairoMakie
include("./Model.jl")

N = 2
C0 = fill(0.01, N)
tspan = (0,2000)

# all initial values
Tr = 0+273.15

Random.seed!(0)
r0 = rand(Normal(0.17, 0.005), N) # from data
α0 =  rand(Normal(-0.28, 1.91), N,N) # from data
α0_diag = rand(truncated(Normal(-0.85, 0.78); upper = 0), N)
α0[diagind(α0)] = α0_diag

# coexistence condition
con1 = r0[2]*α0[1,2] / (r0[1]*α0[2,2])
con2 = r0[2]*α0[1,1] / (r0[1]*α0[2,1])
# dominance condition
con_d = r0[2]*(α0[1,2]+α0[1,1]) / (r0[1]*(α0[2,1]+α0[2,2]))

con1, con2, con_d = log(con1), log(con2), log(con_d)
# T = rand(Uniform(273.15, 25+273.15))
T = 10 + 273.15
ΔT = -1/0.0000862 * (1/T - 1/Tr)

#Base plot
Eα_seq = range(-5, 5, length=100)
Er1 = con1/ΔT .+ Eα_seq
Er2 = con2/ΔT .+ Eα_seq
Erd = con_d/ΔT .+ Eα_seq

colors = ["darkorange" "blue" "black"]
colors_1 = ["darkorange" "blue" "grey"]
plot_s = Plots.plot(Eα_seq, [Er1 Er2 Erd], color = colors, label = ["Survival Condition for i" "Survival Condition for j" "Condition for Dominance"])
Plots.ylims!(-1,1.5)
Plots.xlims!(-5,5)

all_ΔEr = zeros(0)
all_ΔEα = zeros(0)
all_c =  Vector{String}()
for i in 1:200
    Ea_r = rand(Normal(0.95, 0.17), N) # from data
    r = temp_func(T, Tr, r0, Ea_r)#, Ed_r, Th_r) 
    ran_Eα = rand(truncated(Normal(2.19, 3.74); lower = 0, upper = 5), N) # from data
    Ea_α = reshape(repeat(ran_Eα, N),(N,N)) # mean value from data, test variation
    α = temp_func(T, Tr, α0, Ea_α)#, Ed_α, Th_α)

    p = (N = N, α = α, r = r)

    prob = ODEProblem(GLV_model!, C0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))
    # sol.u

    ############# Analysis ##################
    ΔEr = Ea_r[1] - Ea_r[2]
    ΔEα = ran_Eα[1] - ran_Eα[2]
    # ΔE = ΔEr - ΔEα
    div = count(x-> x > 1e-5, sol.u[length(sol)])
    dom = [if div == 1 3 elseif div == 2 && sol.u[length(sol)][1] > sol.u[length(sol)][2] 1 else 2 end]
    append!(all_ΔEr , ΔEr)
    append!(all_ΔEα, ΔEα)
    append!(all_c, colors_1[dom])
end
# Plots.scatter!([all_ΔEα],[all_ΔEr], color = all_c, label = ["i dominates" "j dominates" "cannot coexist"])
Plots.scatter!(all_ΔEα[findall(x->x==colors_1[1], all_c)],all_ΔEr[findall(x->x==colors_1[1], all_c)], color = colors_1[1], label = "i dominates")
Plots.scatter!(all_ΔEα[findall(x->x==colors_1[2], all_c)],all_ΔEr[findall(x->x==colors_1[2], all_c)], color = colors_1[2], label = "j dominates")
Plots.scatter!(all_ΔEα[findall(x->x==colors_1[3], all_c)],all_ΔEr[findall(x->x==colors_1[3], all_c)], color = colors_1[3], label = "cannot coexist")

Plots.plot(plot_s, xaxis = "ΔEα", yaxis = "ΔEr", legend = :topright)
savefig("../../results/Simulations/coex_con.png")

# Plots.plot(sol, xaxis = "Time", yaxis = "Biomass", legend=false)
