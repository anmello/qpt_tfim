import sys
import matplotlib.pyplot as plt
import numpy as np

file1 = str(sys.argv[1])
file2 = str(sys.argv[2])
file3 = str(sys.argv[3])
file4 = str(sys.argv[4])
file5 = str(sys.argv[5])
file6 = str(sys.argv[6])

h1=float(sys.argv[7])
h_val=[-h1,0.,h1]

gs_energies5=np.zeros(3)
gs_magnetization5=np.zeros(3)

gs_energies6=np.zeros(3)
gs_magnetization6=np.zeros(3)

file_one = open(file1, 'r')
file_two = open(file2,'r')
file_three = open(file3,'r')
file_four = open(file4, 'r')
file_five = open(file5, 'r')
file_six = open(file6, 'r')

x1 = np.asarray(file_one.readlines())
y1 = np.asarray(file_two.readlines())
z1 = np.asarray(file_three.readlines())

gs_energies5[0]=x1[0]
gs_energies5[1]=y1[0]
gs_energies5[2]=z1[0]


gs_magnetization5[0]=x1[1]
gs_magnetization5[1]=y1[1]
gs_magnetization5[2]=z1[1]

x2 = np.asarray(file_four.readlines())
y2 = np.asarray(file_five.readlines())
z2 = np.asarray(file_six.readlines())

gs_energies6[0]=x2[0]
gs_energies6[1]=y2[0]
gs_energies6[2]=z2[0]


gs_magnetization6[0]=x2[1]
gs_magnetization6[1]=y2[1]
gs_magnetization6[2]=z2[1]

file_one.close()
file_two.close()
file_three.close()
file_four.close()
file_five.close()
file_six.close()

plt.xlabel('h field')
plt.ylabel('GS energy')
plt.plot(h_val,gs_energies5,label='N=5')
plt.plot(h_val,gs_energies6,label='N=6')
plt.legend()
plt.show()
plt.savefig(fname = 'gs_energy',format='png')

plt.xlabel('h field')
plt.ylabel('GS magnetization')
plt.plot(h_val,gs_magnetization5,label='N=5')
plt.plot(h_val,gs_magnetization6,label='N=6')
plt.legend()
plt.show()
plt.savefig(fname='gs_magnetization',format='png')
