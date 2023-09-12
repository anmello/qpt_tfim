# qpt_tfim
Repository containing codes for the study of a parity dependent quantum phase transition in the transverse field Ising model via VQE.
COMMANDA TO USE
python new_main.py -f 1 -J 1 -B 0 -n 5 -r 2 -i 500

1. "new_main.py" contains the main part of the code, it has to be launched by command prompt specifying flag (0 noiseless simulation, 1 noisy simulation), J (=1 antiferro, =-1 ferro), h (magnetic field) and nqubits (number of qubits in the system). The code outputs a file titled "flag_J_h_nqubits.txt" with the first line containing the GS energy and the second one containing the GS magnetization;
2. "single_plots.py" requires the names of three files (obtained by running "main.py" for each value of h), the number of qubits (nqubits) and a magnetic field value h ( such that {-h,0,h} are the points considered in running "main.py"). The script outputs a plot for the GS energy and a plot for the GS magnetization.
3. "joint_plots.py" requires the names of six files (obtained by running "main.py" for each value of h and for nqubits =5 and nqubits=6) and a magnetic field value h. The script outputs two plots compairing GS energy and GS magnetization in the nqubits=5 and nqubits=6 cases.


