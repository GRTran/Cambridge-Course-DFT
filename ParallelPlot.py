import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('preconditioned_data1.dat')
print(data)

x = [10, 100, 1000, 10000, 100000, 1000000]

legend=['serial','thread=2','thread=3','thread=4','serial vectorised','thread=3 + further parallelisation','thread=3 + vectorisation','thread=2 + further parallelisation wait=active','final optimised']

linestyles1 = ['-','--','-.',':','--','-','--','-.',':']
linemarker1 = ['x','<','.','+','D','>','o','s','*']
#print(data[0:6,0])
for i in range(0,9):
    plt.plot(x,data[6*(i):6*i+6,0],label=legend[i],linestyle=linestyles1[i],marker=linemarker1[i])


plt.legend()
plt.show()
