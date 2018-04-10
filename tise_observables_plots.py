""" tise_obsevables_plots.py

    Creating various plots of observables.

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 4/4/2018
"""
import tise_observables_test as tot
import matplotlib.pyplot as plt
import schrodinger_to_matrix as stm
import numpy as np

## Harmonic Oscillator Basis ##

# Ground State Eigenvalue vs. Increasing No. of Basis Functions #

"""
eigenva_list = []
i_list = []

for i in range(30,1000,100):
    x = tot.discrete_ob_plot([-10,10],i,30,stm.V,511,1,2,0,operator='x')
    eigenva_list.append(x)
    i_list.append(i)

print(eigenva_list)
plt.plot(i_list,eigenva_list)
plt.show()

np.savetxt('2x_grnd_state_10.out',(i_list,eigenva_list))
"""

eigenva_list = []
i_list = []

for i in range(10,50,10):
    x = tot.ho_ob_plot([-1,1],1000,i,stm.V,511,1,1,0,0,operator='x**2')
    eigenva_list.append(x)
    i_list.append(i)
    
print(eigenva_list)
plt.plot(i_list,eigenva_list)

plt.show()

np.savetxt('test1',(i_list,eigenva_list))
