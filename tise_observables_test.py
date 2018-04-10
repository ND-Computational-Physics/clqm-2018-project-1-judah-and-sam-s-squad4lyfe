""" tise_obsevables_test.py

    Finding expectation values of different observables.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/29/2018
"""
import tise_observables_discrete as discrete
import tise_observables_ho as ho
import schrod_plotting as sp
import schrodinger_matrix_solver as sms
import schrodinger_to_matrix as stm 
import numpy as np
import matplotlib.pyplot as plt

def discrete_ob_plot(endpoints,num_steps,dimension,potential,m,hbar,state1,state2,operator='x'):
    dis_test = discrete.Discrete_Observables(endpoints,num_steps,dimension,potential,m,hbar)
    x_val = dis_test.x_set()

    op_matrix = dis_test.gen_matrix_discrete(x_val,operator)
    #print(op_matrix)

    dis_hamilt = sp.schrod_plot_discrete(endpoints, num_steps, potential, 0, 1, m)
    
    eigenstate1 = np.transpose(dis_hamilt[state1])
    eigenstate2 = dis_hamilt[state2]
    
    expectation = np.dot(eigenstate1,np.dot(op_matrix,eigenstate2))
    
    return expectation


def ho_ob_plot(endpoints,num_steps,dimension,potential,m,omega,hbar,state1,state2,operator='x'):
    
    ho_test = ho.HO_Observables(endpoints,num_steps,dimension,potential,m,omega,hbar)
    x_val = ho_test.x_set()

    op_matrix = ho_test.gen_matrix_ho(x_val,operator)
    #print(op_matrix)
    
    ho_hamilt = sp.schrod_plot_ho(endpoints, num_steps, dimension, potential,omega,hbar, 0, 1, m)
    
    eigenstate1 = np.transpose(ho_hamilt[state1])
    eigenstate2 = ho_hamilt[state2]
    
    expectation = np.dot(eigenstate1,np.dot(op_matrix,eigenstate2))
    
    return expectation

    

