# importing datasets
data_r = DataFrame(CSV.File("../../results/TPC/est_params_r.csv"));
data_a = DataFrame(CSV.File("../../results/TPC/est_params_a.csv"));
data_ref = DataFrame(CSV.File("../../results/TPC/ref_biomass.csv"));
data_ref_sd = DataFrame(CSV.File("../../results/TPC/ref_biomass_sd.csv"));
# setting reference temp to 10
r12 = data_r[:,:B0];
Ea_r = data_r[:,:Ea];
r0_all = temp_func((10+273.15), (12+273.15), r12, Ea_r)#, Ed_r, Th_r);

α0_v = -data_a[:,:B0];
row_Eα = data_a[:,:Ea];
Ed_α = fill(4, 3, 3);
row_Th = data_a[:,:Th];
α0_all = temp_func((10+273.15), (12+273.15), α0_v, row_Eα)#, Ed_α, row_Th);

######## W02 + W03 ########
N = 2;
C0 = fill(0.01, N);
tspan = (0,100);
r0 = r0_all[Not(2)];
α0ii = α0_all[Not(2)];
odedata = hcat(data_ref[:,3], data_ref[:,4])
ode_error = hcat(data_ref_sd[:,3], data_ref_sd[:,4]);
time = [1, 24, 49, 72, 91]; # sampled hours 

function LV(du, u, p, t)
    α_ij, α_ji = p
    x, y = u
    du[1] = x*(r0[1]+α0ii[1]*x + α_ij*y)
    du[2] = y*(r0[2]+α0ii[2]*x + α_ji*y)
end

# Take a guess
α_ij = 0.7
α_ji = -0.5
p = [α_ij α_ji];
prob = ODEProblem(LV, C0, tspan, p);
sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=0.5);
# Plots.plot(sol)
# Plots.scatter!(time, odedata)

@model function fitlv(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α_ij ~ truncated(Normal(0.7, 0.5); lower=0, upper=1)
    α_ji ~ truncated(Normal(-0.5, 0.5); lower=-1, upper=0)

    # Simulate Lotka-Volterra model. 
    p = [α_ij,α_ji]
    predicted = solve(prob, AutoTsit5(Rosenbrock23()); p=p, saveat=0.5)

    # Observations.
    if length(predicted) == 201
        for i in 1:length(time)
            t = time[i]
            data[i,:] ~ MvNormal(predicted[t], σ^2 * I)
        end
    end
    return nothing
end
model = fitlv(odedata, prob)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

α_ij = summarystats(chain)[:α_ij,:mean];
α_ji = summarystats(chain)[:α_ji,:mean];
p = [α_ij α_ji];
prob = ODEProblem(LV, C0, tspan, p);
sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=0.5);
α = [α0ii[1] α_ij; α_ji α0ii[2]]

Plots.plot(sol, xaxis = "Time (hours)", yaxis = "Biomass", label="", legend=true, lc = [:darkgreen :darkorange])
Plots.scatter!(time, odedata, yerror=(ode_error, ode_error), label=["S18" "W03"], color = [:darkgreen :darkorange],
 alpha = 0.7, markerstrokecolor=[:darkgreen :darkorange])
# savefig("../../results/Simulations/Tr_S18W03.png")