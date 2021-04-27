from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np

g = 9.81 # m/s^2

def m1_ode(y,t,m1,m2,k1,l01):
    """
    Input:
        - y: Initial conditions, tuple (l1, z1, theta1, z2)
        - t: Time values, array
        - m1,m2,k1,l01: Physical constants, mass of mass 1, mass of mass 2, spring constant of spring 1, equilibrium length of spring 1
    Returns:
        - Derivatives: l1dot, z1dot (l1doubledot), theta1dot, z2dot (theta1doubledot)
    """
    l1, z1, theta1, z2 = y
    l1dot = z1
    theta1dot = z2
    z1dot = 1/2*theta1dot**2 + (m1 + m2)/m1 * g * np.cos(theta1) - k1/m1*(l1 - l01)
    z2dot = -l1dot*z2/l1 - (m1 + m2)/m1 * g * np.sin(theta1)
    return l1dot, z1dot, theta1dot, z2dot

def m2_ode(y,t,m2,k2,l02):
    """
    Input:
        - y: Initial conditions, tuple (l2, z3, theta2, z4)
        - t: Time values, array
        - m2,k2,l02: Physical constants, mass of mass 2, spring constant of spring 2, equilibrium length of spring 2
    Returns:
        - Derivatives: l2dot, z3dot (l2doubledot), theta2dot, z4dot (theta2doubledot)
    """
    l2, z3, theta2, z4 = y
    l2dot = z3
    theta2dot = z4
    z3dot = 1/2*theta2dot**2 + g * np.cos(theta2) - k2/m2*(l2 - l02)
    z4dot = -l2dot*z4/l2 - g * np.sin(theta2)
    return l2dot, z3dot, theta2dot, z4dot

def solve_ode(initial_conditions, t, constants):
    """
    Input:
        - initial_conditions: array [(l10,l1dot0,theta10,theta1dot0),  (l20,l2dot0,theta20,theta2dot0)]
        - t: Time values, array
        - constants: Physical constants, tuple (m1,m2,k1,k2,l01,l02)
    Returns:
        - l1, theta1, l2, theta2: Solutions to all ODEs
    """
    m1,m2,k1,k2,l01,l02 = constants
    spring1_IC, spring2_IC = initial_conditions
    soln1 = odeint(m1_ode, spring1_IC, t, args=(m1,m2,k1,l01))
    soln2 = odeint(m2_ode, spring2_IC, t, args=(m2,k2,l02))
    l1 = soln1[:,0]
    theta1 = soln1[:,2]
    l2 = soln2[:,0]
    theta2 = soln2[:,2]
    return (l1, theta1, l2, theta2)

def convert_to_cartesian(l1, theta1, l2, theta2):
    """
    Input:
        - l1, theta1, l2, theta2: Solution to all ODEs in polar coordinates
    Returns:
        - x1, y1, x1, y2: tuple, Solution to all ODEs in Cartesian coordinates (x1, y1, x2, y2)
    """
    x1 = l1*np.sin(theta1)
    y1 = -l1*np.cos(theta1)
    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)
    return (x1, y1, x2, y2)

print("Pendulum module successfully imported!")