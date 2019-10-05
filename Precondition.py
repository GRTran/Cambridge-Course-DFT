import numpy as np
import matplotlib.pyplot as plt

x = []
x2 = []
y1 = []
y2 = []

for i in range(11,80):
    j = i / 10.0
    x += [j]
    y1 += [1.0/(j-1.0)]

for i in range(0,81):
    j = i/2.0
    x2 += [j]
    y2 += [(8.0 + 4.0*j + 2.0*j**2 + j**3) / (8.0 + 4.0*j + 2.0*j**2 + j**3 + j**4)]

print(x2)

plt.plot(x,y1,label='simple preconditioning function',linestyle='-',marker='o')
plt.plot(x2,y2,label='polynomial preconditioning function',linestyle='--',marker='x')
plt.ylim(0,2.5)
plt.xlim(0,8.0)
plt.legend()
plt.xlabel('Normalised Hamiltonian Diagonal Value')
plt.ylabel('Preconditioning Value')
plt.show()
