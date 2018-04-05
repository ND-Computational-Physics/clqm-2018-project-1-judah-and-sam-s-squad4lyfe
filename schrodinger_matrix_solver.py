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
    """
    Test function to create an (n x n) matrix with random entries.

    Arguments:
            n(int): Dimension of square array
"""

    matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            matrix[i,j] = np.random.random()
    return matrix

def eigensolver(matrix):
    """
    Given a matrix, this function will solve for the eigenvalues and
    eigenfunctions, returning both of these as separate arrays.

    Arguments:
            matrix(array): List of lists representing a matrix
"""

    eigenvalues = np.linalg.eigh(matrix)[0]
    eigenvectors = np.linalg.eigh(matrix)[1]

    return eigenvalues, eigenvectors

if __name__ == "__main__":
    my_matrix = [[1,2,3], [4,5,6], [7,8,9]]
    print(eigensolver(my_matrix)[0])
    print('')
    print(eigensolver(my_matrix)[1])
