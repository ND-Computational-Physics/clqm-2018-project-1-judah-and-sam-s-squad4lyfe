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

def schrod_plot(endpoints, num_points, potential, range_var1, range_var2): 
    """
    Plotting and normalizing our solutions to the Schrodinger Equation.
    
    Inputs:
        endpoints (list): The endpoints at which we want to evaluate between
        num_points (int): The number of points (i.e. solutions) we want to have
        potental (func): The potential we want to consider when evaluating
        range_var1 (int): The lower bound of the eigenvector solutions we want to plot
        range_var2 (int): The upper bound of the eigenvector solutions we want to plot
    
    Outputs:
        (plot): A plot of the normalized eigenvectors you're interested in plotting
    """
    test_case = stm.Schrod_Matrix(endpoints,num_points-1, potential)
    x_val = test_case.x_set()
    matrix = test_case.gen_matrix(x_val)
    
    matrix = .5*matrix # Our matrix is scaled by a factor of 2 too large based on the derivation of this method

    eigenva,eigenve = sms.eigensolver(matrix)

    eigenve = np.transpose(eigenve) # Transposing the eigenvectors so we can plot them against our x_val
    
    h = (test_case.endpoints[1] - test_case.endpoints[0])/test_case.num_steps
    
    for i in range(range_var1,range_var2):
        eigenve_plot = np.array(eigenve[i])*np.sqrt(1/h)
        plt.plot(x_val,eigenve_plot)

    plt.show()

schrod_plot([-5,5],1050,stm.V, 0,2)
