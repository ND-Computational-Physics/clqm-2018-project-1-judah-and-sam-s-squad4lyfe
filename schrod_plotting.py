""" schrod_plotting.py

    Testing and plotting solutions to our Schrodinger solver.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/27/2018
"""
import schrodinger_to_matrix as stm
import numpy as np
import schrodinger_matrix_solver as sms
import matplotlib.pyplot as plt
import scipy
import hermite

def schrod_plot_discrete(endpoints, num_points, potential, range_var1, range_var2, m):
    """
    Plotting and normalizing our solutions to the Schrodinger Equation in the discrete basis.
    
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
    test_case = stm.Schrod_Matrix(endpoints,num_points-1,num_points-1, potential, m)
    x_val = test_case.x_set()
    matrix = test_case.gen_matrix_discrete(x_val)

    eigenva,eigenve = sms.eigensolver(matrix)

    eigenve = np.transpose(eigenve) # Transposing the eigenvectors so we can plot them against our x_val

    h = (test_case.endpoints[1] - test_case.endpoints[0])/test_case.num_steps
    
    for i in range(range_var1,range_var2):
        eigenve_plot = np.array(eigenve[i])*np.sqrt(1/h)
        #plt.plot(x_val,eigenve_plot)
    
    plt.xlabel("(pm)")
    plt.ylabel("(keV)")
    
    #plt.show()
    
    return eigenva, x_val, eigenve_plot

def schrod_plot_ho(endpoints, num_points, dimension, potential, range_var1, range_var2, m):
    """
    Plotting and normalizing our solutions to the Schrodinger Equation in the harmonic oscillator basis.
    
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
    omega = 1 
    
    test_case = stm.Schrod_Matrix(endpoints, num_points-1, dimension, potential, m)
    x_val = test_case.x_set()
    matrix = test_case.gen_matrix_ho(x_val,omega)

    eigenva,eigenve = sms.eigensolver(matrix)

    eigenve = np.transpose(eigenve)

    ho_total = []
    for i in range(0, dimension):
        ho_row = []
        for x in x_val:
            ho_row.append(test_case.ho_soln(x, i, omega, m))
        ho_total.append(ho_row)

    solns = np.matmul(np.transpose(np.array(ho_total)),np.array(eigenve))

    solns = np.transpose(solns)
    
    #for i in range(range_var1, range_var2):
        #plt.plot(x_val, solns[i])
    
    #plt.xlabel("(pm)")
    #plt.ylabel("(keV)")
    #plt.show()
    
    return eigenva, x_val, solns
