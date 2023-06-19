import sys
import matplotlib.pyplot as plt
import numpy as np

file1 = str(sys.argv[1])
file2 = str(sys.argv[2])
file3 = str(sys.argv[3])


file_one = open(file1, 'r')
file_two = open(file2,'r')
file_three = open(file3,'r')

nqubits = int(sys.argv[4])

h1=float(sys.argv[5])
h_val=[-h1,0.,h1]
print(h_val)

gs_energies=np.zeros(3)
gs_magnetization=np.zeros(3)

x = np.asarray(file_one.readlines())
y = np.asarray(file_two.readlines())
z = np.asarray(file_three.readlines())


gs_energies[0]=x[0]
gs_energies[1]=y[0]
gs_energies[2]=z[0]


gs_magnetization[0]=x[1]
gs_magnetization[1]=y[1]
gs_magnetization[2]=z[1]

file_one.close()
file_two.close()
file_three.close()

plt.xlabel('h field')
plt.ylabel('GS energy')
plt.plot(h_val,gs_energies,label='N='+str(nqubits))
plt.legend()
plt.show()
plt.savefig(fname = 'gs_energy',format='png')

plt.xlabel('h field')
plt.ylabel('GS magnetization')
plt.plot(h_val,gs_magnetization,label='N='+str(nqubits))
plt.legend()
plt.show()
plt.savefig(fname='gs_magnetization',format='png')