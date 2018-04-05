""" tise_observables_ho.py

    Finding expectation values of different observables in the harmonic oscillator basis.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/6/2018
"""
import numpy as np
import hermite as hermite
import math
import scipy

class HO_Observables:
    def __init__(self, endpoints, num_steps, dimension, potential, m):
        """
        Arguments:
            endpoints (list): The two endpoints of our range over which we're solving the Schrodinger Equation.

            num_steps (int): The number of steps we're solving the Schrodinger Equation at, essentially the size of our hamiltonian matrix.

            dimension (int): The dimensions of our hamiltonian matrix we'll construct

            potential (function): Our potential, resulting from V(x).

            m (float): The mass of our particle
        """
        self.endpoints = endpoints
        self.num_steps = num_steps
        self.dimension = dimension
        self.potential = potential
        self.m = m
        
    def gen_matrix_ho(self, x, operator='x'):
        """
        Creating our Hamiltonian matrix for the harmonic oscillator basis

        Inputs:
            x (list): A list of our x values we'll evaluate over

        Outputs:
            hamilt (array): Our hamiltonian matrix
        """
        hamilt = np.zeros([self.dimension,self.dimension])

        for i in range(self.dimension):
            for j in range(self.dimension):
                if operator=='x':
                    ele = self.x_element_ho(i,j,x)
                if operator=='x**2':
                    ele = self.x2_element_ho(i,j,x)
                if operator=='p':
                    ele = self.p_element_ho(i,j)
                if operator=='p**2':
                    ele = self.p2_element_ho(i,j)
                hamilt[i][j] = ele

        return hamilt

    def matrix_element_generate_ho(self, i, j, x_points, w):
        """
        Generating the elements of our hamiltonian matrix for a discrete basis.

        Inputs:
            i (index): Row index

            j (index): Column index

            x_points (list): A set of the x points we're evaluating over

            w (float): Our value for omega

        Outputs:
            element (float): The element of our hamiltonian matrix
        """
        h_bar = 1

        prefactor = (-1*h_bar*w)/4

        if i+2 == j:
            element = prefactor*math.sqrt(i+1)*math.sqrt(i+2) + self.ho_integral(x_points, i, j, self.potential, self.m)
        elif i == j:
            element = (-1)*prefactor*((i+1)+(i)) + self.ho_integral(x_points, i, j, self.potential, self.m)
        elif i-2 == j:
            element = prefactor*math.sqrt(i)*math.sqrt(i-1) + self.ho_integral(x_points, i, j, self.potential, self.m)
        else:
            element = 0

        return element

    def ho_integral(self, x_set, i, j, potential, m):
        """
        Evaluating the potential integral in the harmonic oscillator basis.

        Inputs:
            x_set (list): A list of x values we'll evaluate over

            i (int): The row index

            j (int): The column index

            potential (func): Our potential function

            m (float): The mass of our particle

        Outputs:
            returnee (float): The value of our potential integral
        """
        omega = 1
        V = []
        wavefunc_conj = []
        wavefunc = []

        for x in x_set:
            V.append(potential(x))
            wavefunc_conj.append(self.ho_soln(x, i, omega, m))
            wavefunc.append(self.ho_soln(x, j, omega, m))

        solution = np.array(wavefunc_conj)*np.array(V)*np.array(wavefunc)

        returnee = scipy.integrate.simps(solution, x_set)

        return returnee

    def ho_soln(self, x, n, omega, m):
        """
        The solutions to the harmonic oscillator basis, which we'll mulitply by our eigenvectors to the Hamiltonian.

        Inputs:
            x (float): x value we're evaluating at

            n (int): The integer quantum number n

            omega (float): Constant of proprtionality

            m (float): The mass of our particle

        Outputs:
            soln (float): Our solution to part of the HO basis
        """
        h_bar = 1

        x = np.sqrt(m*omega/h_bar)*x
        herm = hermite.hermite(n, x)
        soln = ((m*omega/(np.pi*h_bar))**(1/4))*((math.sqrt((2**n)*math.factorial(n)))**(-1))*herm*np.exp(-.5*x**2)

        return soln

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

# X Operator #

    def V_x(x):
        return x

    def x_element_ho(self,i,j,x_points):
        if i+2 == j:
            element = self.ho_integral(x_points, i, j, V_x, self.m)
        elif i == j:
            element = self.ho_integral(x_points, i, j, V_x, self.m)
        elif i-2 == j:
            element = self.ho_integral(x_points, i, j, V_x, self.m)
        else:
            element = 0

        return element

# X**2 Operator #

    def V_x2(x):
        return x**2

    def x2_element_ho(self,i,j,x_points):
        if i+2 == j:
            element = self.ho_integral(x_points, i, j, HO_Observables.V_x2, self.m)
        elif i == j:
            element = self.ho_integral(x_points, i, j, HO_Observables.V_x2, self.m)
        elif i-2 == j:
            element = self.ho_integral(x_points, i, j, HO_Observables.V_x2, self.m)
        else:
            element = 0

        return element

# P Operator #

    def p_element_ho(self,i,j):
        prefactor = math.sqrt(hbar*m*w/2)

        if i == j+1:
            element = prefator*math.sqrt(j+1)
        elif i == j-1:
            element = -1*prefactor*math.sqrt(j)
        else:
            element = 0
        return element

# P**2 Operator #

    def p2_element_ho(self,i,j):
        
        hbar = 1
        w = 1
        
        prefactor = -1*(hbar*self.m*w/2)

        if i == j+2:
            element = prefactor*math.sqrt((j+1)*(j+2))
        elif i ==j:
            element = prefactor*(-1)*(math.sqrt((j+1)*(j+1))+math.sqrt(j*j))
        elif i == j-2:
            element = prefactor*math.sqrt(j*(j-1))
        else:
            element = 0
        return element
