import numpy as np
import matplotlib.pyplot as plt

precon_dat = np.loadtxt('precondition_presentation_data.dat')
no_precon_dat = np.loadtxt('noprecondition_presentation_data.dat')
print(len(precon_dat))
plt.plot(precon_dat[0:10],precon_dat[10:20], label='preconditioned',linestyle='-',marker='o')
plt.plot(no_precon_dat[0:5],no_precon_dat[5:10],label='no preconditioning',linestyle='--',marker='x')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Plane Waves')
plt.ylabel('Iteration Count')
plt.show()
