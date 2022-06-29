import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pylab as plt
import size_temp_funcs as st
import parameters as par
import model_func as mod


######### Main Code ###########

######## Set up parameters ###########

N = 2 # Number of consumers
M = 2 # Number of resources

# Temperature params
T = 273.15 + 0 # Temperature
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

##### Intergrate system forward #####

def ass_temp_run(N, M, T, Tref, Ma, Ea_D, lf, p_value, typ, K):
    '''
    Main function for the simulation of resource uptake and growth of microbial communities.
    '''
    # Setted Parameters
    k = 0.0000862 # Boltzman constant
    a = 15 # The alpha value for beta distribution in Ea
    B_R0 = 1.70 * np.exp((-0.67/k) * ((1/Tref)-(1/273.15)))/(1 + (0.67/(Ea_D - 0.67)) * np.exp(Ea_D/k * (1/311.15 - 1/Tref))) # Using CUE0 = 0.22, mean growth rate = 0.48
    B_U0 = (1.70/(1 - lf - 0.22)) * np.exp((-0.82/k) * ((1/Tref)-(1/273.15)))/(1 + (0.82/(Ea_D - 0.82)) * np.exp(Ea_D/k * (1/308.15 - 1/Tref)))

    ### Creating empty array for storing data ###
    result_array = np.empty((0,N+M)) # Array to store data in for plotting
    # rich_series = np.empty((0))

    # for i in range(ass):

    ### Resetting values for every assembly ###

    x0 = np.concatenate((np.full([N], 0.1), np.full([M], 1))) # Starting concentration for resources and consumers

    # Set up Ea (activation energy) and B0 (normalisation constant) based on Tom Smith's observations
    B_U = np.random.normal(B_U0, 0.1*B_U0, N) # Adding variation into B0
    B_R = np.random.normal(B_R0, 0.1*B_R0, N)
    T_pk_U = 273.15 + np.random.normal(35, 5, size = N)
    T_pk_R = T_pk_U + 3
    Ea_U = np.random.beta(a, ((a - 1/3) / (0.82/4)) + 2/3 - a, N)*4
    Ea_R = np.random.beta(a, ((a - 1/3) / (0.67/4)) + 2/3 - a, N)*4

    p = np.repeat(p_value, M)  # Resource input

    # Set up model
    U, R, l = par.params(N, M, T, k, Tref, T_pk_U, T_pk_R, B_U, B_R,Ma, Ea_U, Ea_R, np.repeat(Ea_D,N), lf) # Uptake
    l_sum = np.sum(l, axis=1)

    # Integration
    t_fin = 3
    t = np.linspace(0,t_fin-1,t_fin) 
    pars = (U, R, l, p, l_sum, N, M, typ, K) # Parameters to pass onto model

    # pops, infodict = odeint(mod.metabolic_model, y0=x0, t=t, args = pars, full_output=1) # Integrate
    pops = solve_ivp(mod.metabolic_model, t_span= [0,t_fin], y0=x0, t_eval = t, args = pars, method = 'BDF') # Integrate
    pops.y = np.transpose(np.round(pops.y, 7))

    while np.any(np.abs(pops.y[t_fin-1] - pops.y[t_fin-2]) != 0):
        t_fin = t_fin + 1 
        Add_in = solve_ivp(mod.metabolic_model, t_span= [0,t_fin], y0=pops.y[t_fin-2], t_eval = np.linspace(0,1,2), args = pars, method = 'BDF')
        pops.y = np.append(pops.y, [np.transpose(np.round(Add_in.y, 7))[1]], axis = 0)

    # for one resource
    S_e = pops.y[t_fin-1,N:N+M]
    C_e = pops.y[t_fin-1,0:N]
    # alpha_0 = -((1-l_sum)**2*S_e**2)*B_U[0]*B_U[1]/p # when only one resource
    p_i = U[0,0]/np.sum(U, axis = 1)[0]
    p_j = U[1,0]/np.sum(U, axis = 1)[1]
    alpha_0 = -(p_i*p_j*(S_e[0]**2)*(1-lf)*(1-l[0,0])/(1+l[1,0]/(1-l[1,1]))+(1-p_i)*(1-p_j)*(S_e[1]**2)*(1-lf)*(1-l[1,1])/(1+l[0,1]/(1-l[0,0])))*B_U[0]*B_U[1]  # when 2 species and tw oresources
    alpha_E = np.sum(Ea_U)


    alpha = -(U[0,0]*U[1,0]*(1-lf)*S_e[0]/(C_e[0]*U[0,0]+C_e[1]*U[1,0])+U[0,1]*U[1,1]*(1-lf)*S_e[1]/(C_e[0]*U[0,1]+C_e[1]*U[1,1]))
    ### Storing simulation results ###
    result_array = np.append(result_array, pops.y, axis=0)

    return result_array, alpha_0, alpha_E, alpha


# result_array, alpha_0, alpha_E, alpha = ass_temp_run(N, M, T, Tref, Ma, Ea_D, lf, p_value, typ, K)

# print(alpha)
# print(alpha_0)
