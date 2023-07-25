#!/bin/bash
#   f 0 noiseless simulation, 1 noisy simulation, 2 real backend), J (=1 antiferro, =-1 ferro), B (magnetic field) and n (number of qubits in the system)
f=1
J=1
n=6
r=2
i=500

for B in -0.3 0 0.3
do
	python new_main.py -f $f -J $J -B $B -n $n -r $r -i $i > output_log/output_${f}_${B}_${J}_${n}_${r}_${i}.txt &
done
