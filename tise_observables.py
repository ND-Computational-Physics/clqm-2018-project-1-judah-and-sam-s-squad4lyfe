""" tise_observables.py

    Finding expectation values of different observables.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/6/2018
"""

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
