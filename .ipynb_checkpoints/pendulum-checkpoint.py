from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

g = 9.8 # m/s^2

def x1y1_odes(y, t, k1, m1, k2, m2, l01, l02 ):
    x1, y1, z1, z2 = y
    x1dot = z1
    z1dot = (k1*x1/m1) * (1 - l01/np.sqrt(x1**2 + y1**2))
    y1dot = z2
    z2dot = (m1+m2)*g/m2 + (k1*y1/m2) * (1 - l01/np.sqrt(x1**2 + y1**2))
    return x1dot, y1dot, z1dot, z2dot

def x2y2_odes(y,t, k1, m1, k2, m2, l01, l02 ):
    x2, y2, z3, z4 = y
    x2dot = z3
    z3dot = (k2*x2/m2) * (1 - l02/np.sqrt(x2**2 + y2**2))
    y2dot = z4
    z4dot = g + (k2*y2/m2) * (1 - l02/np.sqrt(x2**2 + y2**2))
    return x2dot, y2dot, z3dot, z4dot

def solve_odes(spring1_IC, spring2_IC, t, spring_constants):
    k1, m1, k2, m2, l01, l02 = spring_constants 
    x1y1_soln = odeint(x1y1_odes, spring1_IC, t, args=(k1, m1, k2, m2, l01, l02 ))
    print(x1y1_soln)
    x2y2_soln = odeint(x2y2_odes, spring2_IC, t, args=(k1, m1, k2, m2, l01, l02 ))
    return (x1y1_soln, x2y2_soln)

def plot_soln(x1, y1, x2, y2):
    plt.plot(x1, y1)
    plt.plot(x2, y2)

print("Pendulum module successfully imported!")

