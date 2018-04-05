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

def discrete_ob_plot(endpoints,num_steps,dimension,potential,m,range_var1,range_var2,operator='x'):
    dis_test = discrete.Discrete_Observables(endpoints,num_steps,dimension,potential,m)
    x_val = dis_test.x_set()

    op_matrix = dis_test.gen_matrix_discrete(x_val,operator)
    #print(op_matrix)

    eigenva,eigenve = sms.eigensolver(op_matrix)

    eigenve = np.transpose(eigenve) # Transposing the eigenvectors so we can plot them against our x_val

    h = (dis_test.endpoints[1] - dis_test.endpoints[0])/dis_test.num_steps
    
    for i in range(range_var1,range_var2):
        eigenve_plot = np.array(eigenve[i])*np.sqrt(1/h)
        plt.plot(x_val,eigenve_plot)
    
    plt.xlabel("(pm)")
    plt.ylabel("(keV)")
    
    plt.show()
    
    return eigenva, x_val, eigenve_plot

#discrete_ob_plot([-1,1],100,30,stm.V,511,0,10,operator='x')

def ho_ob_plot(endpoints,num_steps,dimension,potential,m,omega,hbar,range_var1,range_var2,operator='x'):
    
    ho_test = ho.HO_Observables(endpoints,num_steps,dimension,potential,m,omega,hbar)
    x_val = ho_test.x_set()

    op_matrix = ho_test.gen_matrix_ho(x_val,operator)
    #print(op_matrix)
    
    eigenva,eigenve = sms.eigensolver(op_matrix)

    eigenve = np.transpose(eigenve)

    ho_total = []
    for i in range(0, dimension):
        ho_row = []
        for x in x_val:
            ho_row.append(ho_test.ho_soln(x, i, omega, m))
        ho_total.append(ho_row)

    solns = np.matmul(np.transpose(np.array(ho_total)),np.array(eigenve))

    solns = np.transpose(solns)
    
    for i in range(range_var1, range_var2):
        plt.plot(x_val, solns[i])
    
    plt.xlabel("(pm)")
    plt.ylabel("(keV)")
    plt.show()
    
    return eigenva

ans = ho_ob_plot([-1,1],1000,30,stm.V,511,1,1,0,3,operator='p**2')
print(ans)
