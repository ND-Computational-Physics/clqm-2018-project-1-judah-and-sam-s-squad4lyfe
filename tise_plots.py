""" tise_plots.py

    Actual plots of analysis of the results from our time-independent Schrodinger equation solver

    Language: Python 3

    Sam Porter, Judah Van Zandt
    University of Notre Dame

    Written for Computational Lab in Quantum Mechanics, Spring 2018
    Last updated on 3/27/2018
"""

#x = schrod_plot_discrete([-1,1],1000,stm.V,0,6,511)
#y = schrod_plot_ho([-1,1],1000,30,stm.V,0,6,511)

"""
# Discrete Energy vs. Quantum Number Plot while Changing Number of Solutions #


for i in range(10,1000,100):
    x = schrod_plot_discrete([-1,1],i,stm.V,0,1,511)

    plt.plot(range(0,i-2),x[0])
    plt.xlabel("n")
    plt.ylabel("(keV)")

plt.show()


# Discrete Energy vs. Quantum Number Plot while Changing Range #

for i in range(1,10):
    x = schrod_plot_discrete([-i,i],1000,stm.V,0,1,511)

    plt.plot(range(0,998),x[0])
    plt.xlabel("n")
    plt.ylabel("(keV)")

plt.show()


# Discrete Ground State Energies vs. Diff. Dimensions #

grnd_energy = []
for step in range(10, 200, 10):
    x = schrod_plot_discrete([-1,1],step,stm.V,0,1,511)
    eigenva = x[0]
    grnd_energy.append(eigenva[0])
x = list(range(10,200,10))
y = grnd_energy

plt.plot(x,y)
plt.xlabel("Number of Eigenfunctions")
plt.ylabel("(keV)")
plt.show()
"""

# Harmonic Oscillator Ground State Energies vs. Diff. Number of Steps #

grnd_energy = []
for step in range(10, 200, 10):
    x = schrod_plot_ho([-1,1],step,30,stm.V,0,1,511)
    eigenva = x[0]
    grnd_energy.append(eigenva[0])
x = list(range(10,200,10))
y = grnd_energy

print(grnd_energy)
plt.plot(x,y)
plt.xlabel("Number of Steps")
plt.ylabel("(keV)")
plt.show()

"""
# Harmonic Oscillator Energy vs. Quantum Number Plot #

x = schrod_plot_ho([-1,1],1000,30,stm.V,0,1,511)

plt.plot(range(0,30),x[0])
plt.xlabel("n")
plt.ylabel("(keV)")
plt.show()

"""
