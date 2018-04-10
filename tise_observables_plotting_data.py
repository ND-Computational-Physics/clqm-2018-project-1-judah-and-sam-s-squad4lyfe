""" tise_obsevables_plotting_data.py

    Creating various plots of observables, taking data we had saved from .out files to do so.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 4/4/2018
"""
import numpy as np
import matplotlib.pyplot as plt

for w in ['.25',1]:
    (i,eigenva) = np.loadtxt('test{0}'.format(str(w)))

    plt.plot(i,np.abs(eigenva),label="w = {0}".format(str(w)))

plt.xlabel('Dimension')
plt.ylabel('<2|x**2|0>')
plt.legend(loc="upper right")
plt.show()
