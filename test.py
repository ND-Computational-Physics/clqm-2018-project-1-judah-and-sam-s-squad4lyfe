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

test_case = stm.Schrod_Matrix([0,100],10,stm.V)

x_val = test_case.x_set()
matrix = test_case.gen_matrix(x_val)
matrix = .5*matrix

print(matrix)
