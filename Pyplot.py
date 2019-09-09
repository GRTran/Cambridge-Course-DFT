import numpy as np
import matplotlib.pyplot as plt

wavefunction = np.loadtxt('wvfn_1.dat')
localpotential = np.loadtxt('pot.dat')

plt.plot(wavefunction[:,0], wavefunction[:,1])
plt.show()
plt.plot(localpotential[:,0], localpotential[:,1])
plt.show()
