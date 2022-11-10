import numpy as np
from scipy.integrate import solve_ivp
r = np.array((0.1, 0.4))
a = np.array(([-0.2,-0.1],[-0.1,-0.2]))
pars = (r, a)
t_fin = 100
t = np.linspace(0,t_fin-1,t_fin) 
x0 = np.full(2, 0.1)

def lv(t, x, r, a):
    dx = x * (r+a@x)
    return dx


pops = solve_ivp(lv, t_span= [0,t_fin],
 y0=x0, t_eval = t, args = pars, method = 'BDF')

np.transpose(pops.y)[t_fin-1]