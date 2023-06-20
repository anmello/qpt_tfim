#!/bin/bash
#   f 0 noiseless simulation, 1 noisy simulation, 2 real backend), J (=1 antiferro, =-1 ferro), B (magnetic field) and n (number of qubits in the system)
J=-1
for B in -0.2 0 0.2
do
	python main.py -f 0 -J $J -B $B -n 5  > output_log/output_${B}_${J}.txt &
done
