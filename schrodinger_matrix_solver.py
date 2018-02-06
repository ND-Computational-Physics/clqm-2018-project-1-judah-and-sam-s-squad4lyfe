""" schrodinger_matrix_solver.py
Language: Python3
Sam Porter, Judah Van Zandt
January 28, 2018
Written for CLQM Schrodinger Project
"""

# The basis of this project is to convert the linear schrodinger eqaution into
# an eigenvalue problem. This means we have to set up a modifiable Hamiltonian
# matrix, as well as a way to solve for its eigenvalues.
import numpy as np
np.set_printoptions(precision=3)

def random_matrix(n):

    matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            matrix[i,j] = np.random.random()
    return matrix

def eigensolver(matrix):

    eigenvalues = np.linalg.eig(matrix)[0]
    eigenvectors = np.linalg.eig(matrix)[1]

    return eigenvalues, eigenvectors

if __name__ == "__main__":
    print(eigensolver(random_matrix(3))[0])
    print('')
    print(eigensolver(random_matrix(3))[1])
