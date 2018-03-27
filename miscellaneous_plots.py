"""
Miscellaneous Plotting Functions
Language: python3
Judah Van Zandt
Computational Lab in Quantum Mechanics
March 9, 2018
"""

import schrodinger_to_matrix as stm
import numpy as np
import schrodinger_matrix_solver as sms
import matplotlib.pyplot as plt
import scipy
import hermite
import schrod_plotting

def infinite_well_potential(x):
    potential = 0
    return potential

def harmonic_oscillator_potential(x):
    potential = .5*511*x**2
    return potential


def discrete_energy_vs_n(title,lower_bound, upper_bound, step_size=1, V = stm.V):
    """
    This function plots the energies of the excited states E_n as a function of E.
    For the HO potential, we expect this plot to be linear, while for the infinite
    potential well, we expect it to be quadratic.
    """
    for i in range(lower_bound,upper_bound,step_size):
        x=range(i-2)
        plt.plot(x,schrod_plotting.schrod_plot_discrete([-1,1],i,V,0,1,511))
    plt.xlabel("n")
    plt.ylabel("Energy of the nth Excited State")
    plt.title("Energies of the {} vs. n ".format(title))
#discrete_energy_vs_n("Infinite Square Well", 999, 1000, 1, infinite_well_potential)

def harmonic_energy_vs_n(endpoints, numpoints, dimension,potential, range_var1, range_var2, m=511):
    """
    This function plots the energies of the excited states E_n as a function of E.
    For the HO potential, we expect this plot to be linear, while for the infinite
    potential well, we expect it to be quadratic.
    """
    eigenva = schrod_plotting.schrod_plot_ho(endpoints, numpoints, dimension,potential, range_var1, range_var2, 511)
    x_values = np.linspace(0, len(eigenva)-1, len(eigenva))
    y_values = eigenva

    plt.plot(x_values, y_values)
    plt.xlabel("n")
    plt.ylabel("Energy of the nth Excited State")
    plt.title("Energies of the Harmonic Oscillator vs. n")


def discrete_energy_ground(max_eig, V = stm.V):
    """
    This function is for making a plot of the value of the ground state energy
    vs. the number of basis vecotr used to calculate E. So as we move along the
    x-axis, we expect that the function will converge to the proper value of E_0.
    """
    x_values = []
    y_values = []
    for i in range(3,max_eig):
        x_values.append(i)
        y_values.append(schrod_plotting.schrod_plot_discrete([-1,1],i,V,0,1,511)[0])
    plt.plot(x_values, y_values)
    plt.xlabel("Number of Basis Functions")
    plt.ylabel("Ground State Energy")
    print(x_values)
    print(y_values)

def harmonic_energy_ground(max_eig, V = stm.V):
    """
    This function is for making a plot of the value of the ground state energy
    vs. the number of basis vecotr used to calculate E. So as we move along the
    x-axis, we expect that the function will converge to the proper value of E_0.
    """
    x_values = []
    y_values = []
    for i in range(3,max_eig):
        x_values.append(i)
        y_values.append(schrod_plotting.schrod_plot_ho([-1,1],100,i,stm.V,0,1,511)[0])
    plt.plot(x_values, y_values)
    plt.xlabel("Number of Basis Functions")
    plt.ylabel("Ground State Energy")
    print(x_values)
    print(y_values)



#harmonic_energy_vs_n([-1,1], 100, 5,stm.V, 0, 1, m=511)
discrete_energy_ground(100, V = stm.V)

plt.show()
