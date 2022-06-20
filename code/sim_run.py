from Community_model import ass_temp_run
import matplotlib.pylab as plt
import numpy as np

N = 2 # Number of consumers
M = 1 # Number of resources

# Temperature params
Tref = 273.15 + 0 # Reference temperature Kelvin
Ma = 1 # Mass
Ea_D = 3.5 # Deactivation energy - only used if use Sharpe-Schoolfield temp-dependance
lf = 0.4 # Leakage
p_value = 1 # External input resource concentration

# Assembly
# ass = 1 # Assembly number, i.e. how many times the system can assemble
# t_fin = 2000 # Number of time steps
typ = 1 # Functional response, Type I or II
K = 5 # Half saturation constant
# strc = 1
k = 0.0000862 # Boltzman constant

alpha_all = np.empty((0))
##### Intergrate system forward #####
for i in range(6):
    T = 273.15 + 5*i # Temperature
    np.random.seed(0)
    result_array, alpha_0, alpha_E, alpha = ass_temp_run(N, M, T, Tref, Ma, Ea_D, lf, p_value, typ, K)
    alpha_all = np.append(alpha_all, alpha)

T_plot = np.array(range(0, 30, 5))
y = alpha_0[0] * np.exp((-alpha_E/k) * ((1/T_plot)-(1/Tref)))