""" schrod_plotting.py

    Testing and plotting solutions to our Schrodinger solver.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 2/8/2018
"""
import schrodinger_to_matrix as stm
import numpy as np
import schrodinger_matrix_solver as sms
import matplotlib.pyplot as plt
import scipy
import hermite

def schrod_plot_discrete(endpoints, num_points, potential, range_var1, range_var2, m): 
    """
    Plotting and normalizing our solutions to the Schrodinger Equation.
    
    Inputs:
        endpoints (list): The endpoints at which we want to evaluate between
        num_points (int): The number of points (i.e. solutions) we want to have
        potental (func): The potential we want to consider when evaluating
        range_var1 (int): The lower bound of the eigenvector solutions we want to plot
        range_var2 (int): The upper bound of the eigenvector solutions we want to plot
        m (float): The mass of our particle
    
    Outputs:
        (plot): A plot of the normalized eigenvectors you're interested in plotting
    """
    test_case = stm.Schrod_Matrix(endpoints,num_points-1, potential, m)
    x_val = test_case.x_set()
    matrix = test_case.gen_matrix(x_val)

    eigenva,eigenve = sms.eigensolver(matrix)

    eigenve = np.transpose(eigenve) # Transposing the eigenvectors so we can plot them against our x_val
    
    h = (test_case.endpoints[1] - test_case.endpoints[0])/test_case.num_steps
    
    for i in range(range_var1,range_var2):
        eigenve_plot = np.array(eigenve[i])*np.sqrt(1/h)
        plt.plot(x_val,eigenve_plot)
    
    plt.xlabel("(pm)")
    plt.ylabel("(keV)")
    
    plt.show()
    return eigenva

def ho_soln(x, n, omega, m):
    """
    The solutions to the harmonic oscillator basis, which we'll mulitply by our eigenvectors to the Hamiltonian.
    
    Inputs:
        x (float): x value we're evaluating at
        n (int): The integer quantum number n
        omega (float): Constant of proprtionality
        m (float): The mass of our particle
    Outputs:
        soln (float): Our solution to part of the HO basis
    """
    x = np.sqrt(m*omega)*x
    herm = hermite.hermite(n, x)
    soln = ((m*omega/np.pi)**(1/4))(1/np.sqrt((2**n)*np.factorial(n)))*herm*np.exp(-.5*x**2)
    
    return soln

def schrod_plot_ho(endpoints, num_points, potential, range_var1, range_var2, m):
    
    eigenva,eigenve = sms.eigensolver(matrix)

    eigenve = np.transpose(eigenve)
    
    omega = 1 
    
    ho = []
    for i in range(range_var1, range_var2):
        for x in x_val:
            ho.append(ho_soln(x, i, omega, m))
        eigenve_plot = np.array(ho)*np.array(eigenve[i])
        plt.plot(x_val,eigenve_plot)
    
    plt.xlabel("(pm)")
    plt.ylabel("(keV)")
    plt.show()

x = schrod_plot_discrete([-2,2],1000,stm.V, 0,3, 511)
print(x)
