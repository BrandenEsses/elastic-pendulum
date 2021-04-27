from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

g = 9.81 # m/s^2

def m1_ode(y,t,m1,m2,k1,k2,l01,l02):
    l1, z1, theta1, z2 = y
    l1dot = z1
    theta1dot = z2
    z1dot = 1/2*theta1dot**2 + (m1 + m2)/m1 * g * np.cos(theta1) - k1/m1*(l1 - l01)
    z2dot = -l1dot*z2/l1 - (m1 + m2)/m1 * g * np.sin(theta1)
    return l1dot, z1dot, theta1dot, z2dot

def m2_ode(y,t,m1,m2,k1,k2,l01,l02):
    l2, z3, theta2, z4 = y
    l2dot = z3
    theta2dot = z4
    z3dot = 1/2*theta2dot**2 + g * np.cos(theta2) - k2/m2*(l2 - l02)
    z4dot = -l2dot*z4/l2 - g * np.sin(theta2)
    return l2dot, z3dot, theta2dot, z4dot

def solve_ode(initial_conditions, t, constants):
    m1,m2,k1,k2,l01,l02 = constants
    spring1_IC, spring2_IC = initial_conditions
    soln1 = odeint(m1_ode, spring1_IC, t, args=(m1,m2,k1,k2,l01,l02))
    soln2 = odeint(m2_ode, spring2_IC, t, args=(m1,m2,k1,k2,l01,l02))
    l1 = soln1[:,0]
    theta1 = soln1[:,2]
    l2 = soln2[:,0]
    theta2 = soln2[:,2]
    return (l1, theta1, l2, theta2)

print("Pendulum module successfully imported!")

