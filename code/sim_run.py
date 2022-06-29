# from cProfile import label
from Community_model import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

N = 2 # Number of consumers
M = 2 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant
k = 0.0000862 # Boltzman constant

for rs in range(3):

    alpha_all = np.empty((0))
    alpha_all_0 = np.empty((0))
    results_all = np.empty((0,2))
    ##### Intergrate system forward #####

    for i in range(5):
        T = 273.15 + 5*i # Temperature
        np.random.seed(rs)
        result_array, alpha_0, alpha_E, alpha = ass_temp_run(N, M, T, Tref, Ma, Ea_D, lf, p_value, typ, K)
        alpha_all_0 = np.append(alpha_all_0, alpha_0)
        alpha_all = np.append(alpha_all, alpha)
        results_all = np.append(results_all, [result_array[len(result_array)-1,0:N]], axis = 0)

    T_plot = np.array(range(0, 25, 5))+273.15
    y = np.mean(alpha_all_0) * np.exp((-alpha_E/k) * ((1/T_plot)-(1/Tref)))
    
    print("Equilibrium biomass at all temps:\n", results_all)
    print("α0 at all temps:" , alpha_all_0, "; mean:", np.mean(alpha_all_0))
    print("α at all temps:" , alpha_all)
    print("Estimation by Arrhenius:", y)

    plt.plot(T_plot - 273.15, alpha_all, label = "alpha")
    plt.plot(T_plot - 273.15, y, label = "estimation")
    plt.xlabel("T")
    plt.title(f"Scenario %d" %(rs+1))
    plt.legend()
    plt.show()