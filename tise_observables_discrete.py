""" tise_observables_discrete.py

    Finding expectation values of different observables in the discrete basis.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/6/2018
"""
class Discrete_Observables:
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
        
    def gen_matrix_discrete(self, x):
        """
        Creating our Hamiltonian matrix for the discrete basis

        Inputs:
            x (list): A list of our x values we'll evaluate over

        Outputs:
            hamilt (array): Our hamiltonian matrix
        """
        hamilt = np.zeros([self.num_steps-1,self.num_steps-1])

        for i in range(self.num_steps-1):
            for j in range(self.num_steps-1):
                ele = self.matrix_element_generate_discrete(i,j,x)
                hamilt[i][j] = ele

        return hamilt
    
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
    
    def x_element_discrete(self,i,j,x_points):
        if i == j:
            element = x_points[i]
        else:
            element = 0
        
        return element

# X**2 Operator #

    def x2_element_discrete(self,i,j,x_points):
        if i == j:
            element = (x_points[i])**2
        else:
            element = 0
    
    return element

# P Operator #

# P**2 Operator #

    def p2_element_discrete(self,i,j,x_points):
        hbar = 1

        h = (self.endpoints[1] - self.endpoints[0])/self.num_steps

        if i == j:
            element = (hbar**2)*(2/h**2)
        elif i >= j+2:
            element = 0
        elif j >= i+2:
            element = 0
        else:
            element = -1/(h**2)

        return element
