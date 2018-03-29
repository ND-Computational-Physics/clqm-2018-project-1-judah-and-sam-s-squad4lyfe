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

dis_test = discrete.Discrete_Observables([-1,1],100,30,10,511)
x_val = dis_test.x_set()

op_matrix = dis_test.gen_matrix_discrete(x_val,operator='p**2')
#print(op_matrix)

ho_test = ho.HO_Observables([-1,1],100,30,10,511,1)
x_val = ho_test.x_set()

op_matrix = ho_test.gen_matrix_ho(x_val,operator='x')
print(op_matrix)
