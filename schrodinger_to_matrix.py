""" schrodinger_to_matrix.py

    Rewriting the Schrodinger Equation as a matrix diagonalization problem.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 2/8/2018
"""
import numpy as np
import schrod_plotting as sp

def V(x):
    """
    Our potential function we want to find solution to the Schrodinger Equation for.

    Arguments:
        x (float): Discrete values between the two endpoints of evaluation.
    """
    V_func = 0 # This is the potential function! Change it here!

    return V_func


class Schrod_Matrix:
    def __init__(self, endpoints, num_steps, potential, m):
        """
        Arguments:
            endpoints (list): The two endpoints of our range over which we're solving the Schrodinger Equation.

            num_steps (int): The number of steps we're solving the Schrodinger Equation at, essentially the size of our hamiltonian matrix.

            potential (function): Our potential, resulting from V(x).
            
            m (float): The mass of our particle
        """
        self.endpoints = endpoints
        self.num_steps = num_steps
        self.potential = potential
        self.m = m

    def matrix_element_generate(self, i, j, x_points):
        """
        Generating the elements of our hamiltonian matrix.

        Arguments:
            i (index): Row index

            j (index): Column index

            x_points (list): A set of the x points we're evaluating over
        """
        hbar = 1
        
        h = (self.endpoints[1] - self.endpoints[0])/self.num_steps
        
        if i == j:
            element = ((hbar**2)/(2*self.m))*(2/h**2) + self.potential(x_points[i])
        elif i >= j+2:
            element = 0
        elif j >= i+2:
            element = 0
        else:
            element = (1/(2*self.m))*(-1/(h**2))

        return element

    def x_set(self):
        """
        Defining our set of x values which we want to evaluate over.
        
        Outputs:
            x (list): A list of our x values
        """
        x = []
        step = (self.endpoints[1] - self.endpoints[0])/self.num_steps

        for i in range(1,self.num_steps):
            x.append(self.endpoints[0] + i*step)
        
        return x

    def gen_matrix(self,x):
        """
        Creating our Hamiltonian matrix
        
        Inputs:
            x (list): A list of our x values we'll evaluate over
        
        Outputs:
            hamilt (array): Our hamiltonian matrix
        """
        hamilt = np.zeros([self.num_steps-1,self.num_steps-1])

        for i in range(self.num_steps-1):
            for j in range(self.num_steps-1):
                ele = self.matrix_element_generate(i,j,x)
                hamilt[i][j] = ele

        return hamilt
    
    def ho_integral(x_set, i, j, potential, m):
        omega = 1
        V = []
        wavefunc_conj = []
        wavefunc = []
        
        for x in x_set:
            V.append(potential(x))
            wavefunc_conj.append(sp.ho_soln(x, i, omega, m))
            wavefunc.append(sp.ho_soln(x, j, omega, m))
        
        solution = np.array(wavefunc_conj)*np.array(V)*np.array(wavefunc)
        
        return np.trapz(solution, x_set)
            
