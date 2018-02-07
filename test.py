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



x_val = test_case.x_set()
matrix = test_case.gen_matrix(x_val)
matrix = .5*matrix

eigenva,eigenve = sms.eigensolver(matrix)

print(matrix)
print(eigenva)
print(eigenve)


#plt.plot(x_val,stm.V(np.array(x_val)))

for i in range(len(eigenva)):
    plt.plot(x_val,eigenve[i])

plt.show()
