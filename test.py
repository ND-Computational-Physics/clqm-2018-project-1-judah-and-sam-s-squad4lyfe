""" test.py

    Testing our Schrodinger solver.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 2/1/2018
"""
import schrodinger_to_matrix as stm
import numpy as np
import schrodinger_matrix_solver as sms
import matplotlib.pyplot as plt

test_case = stm.Schrod_Matrix([-5,5],200, stm.V)
x_val = test_case.x_set()
matrix = test_case.gen_matrix(x_val)
matrix = .5*matrix

eigenva,eigenve = sms.eigensolver(matrix)

eigenve = np.transpose(eigenve)

#plt.plot(x_val,stm.V(np.array(x_val)))
x=[]
y=[]
print(len(eigenve))
print(len(x_val))
for i in range(len(eigenve)):
    x.append(x_val[i])
    y.append(eigenve[0][i])
print(x)
print(y)
plt.plot(x,y)

for i in range():
    eigenve1 = np.array(eigenve[i])**2
    norm_factor = np.trapz(eigenve1,x_val)
    eigenve_plot = np.array(eigenve[i])*np.sqrt(1/norm_factor)
    plt.plot(x_val,eigenve_plot)

plt.show()
