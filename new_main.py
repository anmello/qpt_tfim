import numpy as np 
from qiskit_ibm_runtime import QiskitRuntimeService, Options, Session, Estimator
from qiskit.circuit.library import EfficientSU2, RealAmplitudes
from qiskit.algorithms.optimizers import SPSA
from qiskit.algorithms.minimum_eigensolvers import VQE

from qiskit.quantum_info import Pauli, SparsePauliOp

#hamiltonian generation (for N=5,6)
# number of qubits = i+1+(N-i-1)=N qubits
# see https://qiskit.org/documentation/migration_guides/opflow_migration.html#qubit-paulis
from qiskit.quantum_info import Pauli, SparsePauliOp
from qiskit.algorithms.minimum_eigensolvers import NumPyMinimumEigensolver
from qiskit.providers.fake_provider import FakeAuckland
from qiskit_aer.noise import NoiseModel

import argparse
import datetime
'''
   Read from config file.
   USAGE example: python main.py -f 0 -J 1 -B 0.1 -n 5
   f 0 noiseless simulation, 1 noisy simulation, 2 real backend), J (=1 antiferro, =-1 ferro), h (magnetic field) and nqubits (number of qubits in the system)
 '''

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--flagnoise', type=int, required=True)
parser.add_argument('-J', '--J', type=float, required=True)
parser.add_argument('-B', '--magneticfield', type=float, required=True)
parser.add_argument('-n', '--nqubits', type=int, required=True)
parser.add_argument('-r', '--reps', type = int, required= True)
parser.add_argument('-i', '--maxiter', type = int, required= True)

args = parser.parse_args()

print(args.flagnoise, args.J, args.magneticfield, args.nqubits)

flag = int(args.flagnoise)

J = float(args.J)

h = float(args.magneticfield)

nqubits = int(args.nqubits)

nreps = int(args.reps)

niter = int(args.maxiter)

def ham_generator(N, h, J):
    l = [] #we'll add the terms (term,coef)
    # field term
    for i in range(N):
        t = "I"*i + "Z" + "I"*(N-1-i)
        l.append( (t,h) )
    # interaction term
    for k in range(N - 1):
        t = "I"*k + "XX" + "I"*(N-2-k)
        l.append( (t,J) )
    # PBC
    t = "X" + "I"*(N-2) + "X"
    l.append( (t,J) )
    H = SparsePauliOp.from_list(l)
    return H
    
# definition magnetization
def magnetization(N):
    l=[]
    coef=1/N
    # field term
    for i in range(N):
        t = "I"*i + "X" + "I"*(N-1-i)
        l.append( (t,coef) )
    M = SparsePauliOp.from_list(l)
    return M

#DEFINITION OF HAMILTONIAN AND OPERATOR

H = ham_generator(nqubits, h, J)
M = magnetization(nqubits)


sol = NumPyMinimumEigensolver().compute_minimum_eigenvalue(H, aux_operators=[M])

print('correct solution is: ', sol.eigenvalue)

service = QiskitRuntimeService(channel='ibm_quantum')

sim = "ibmq_qasm_simulator"

ansatz = RealAmplitudes(num_qubits=nqubits, reps=nreps)
#ansatz = EfficientSU2(nqubits, reps=3, entanglement='linear', insert_barriers=False)
#ansatz_true = EfficientSU2(nqubits, reps=3, entanglement='linear', insert_barriers=True)

optimizer = SPSA(maxiter=niter)

np.random.seed(6)
initial_point = np.random.random(ansatz.num_parameters)

# Create an object to store intermediate results
from dataclasses import dataclass

@dataclass
class VQELog:
    values: list
    parameters: list

    def update(self, count, parameters, mean, _metadata):
        self.values.append(mean)
        self.parameters.append(parameters)
        print(f"Running circuit {count} of ~400", end="\r", flush=True)

log = VQELog([], [])

backend = service.backend("ibm_cairo")
# Set options to include the noise model
options = Options()

# Set number of shots, optimization_level and resilience_level
options.execution.shots = 1024
options.optimization_level = 3
options.resilience_level = 1


if flag == 0:
    print('NOISELESS')
    options.simulator = {
        "seed_simulator": 42
    }
    with Session(service=service, backend=sim) as session:
        estimator = Estimator(session=session, options=options)
        vqe = VQE(estimator=estimator,
                  ansatz=ansatz,
                  optimizer=optimizer,
                  callback=log.update,
                  initial_point=initial_point)
        M = magnetization(nqubits)
        H = ham_generator(nqubits, h, J)
        result = vqe.compute_minimum_eigenvalue(H, aux_operators=[M])

        print("Experiment complete.".ljust(30))
        print(f"GS: {result.optimal_value}")
        print(f"Magn:{result.aux_operators_evaluated}")
        session.close()

    print('noisy result for the gs is: ', result.eigenvalue)

    tmp_title = "flag" + "_" + str(J) + "_" + str(h) + "_"  + "_" +str(nqubits)+ "_" + str(nreps) + "_" +str(niter)+".txt"
    values = open(tmp_title, "w+")
    print(result.optimal_value, file=values)
    print((result.aux_operators_evaluated)[0][0], file=values)

    values.close()

    conv_title = str(flag)+"_"+str(J)+"_"+str(h) + "_" + "_" + str(nqubits) + "_"  + str(nreps) + "_"+str(niter) + "_"+ "conv" + ".txt"
	
    conv = open(conv_title,"w+")
    print(log.values, file = conv)
	
    conv.close()

elif flag == 1:
    print('NOISY')
    noise_model = NoiseModel.from_backend(backend)
    basis_gates = backend.configuration().basis_gates
    coupling_map = backend.configuration().coupling_map

    options.simulator = {
        "noise_model": noise_model,
        "basis_gates": basis_gates,
        "coupling_map": coupling_map,
        "seed_simulator": 42
    }
    with Session(service=service, backend=sim) as session:
        estimator = Estimator(session=session, options=options)
        vqe = VQE(estimator=estimator,
                  ansatz=ansatz,
                  optimizer=optimizer,
                  callback=log.update,
                  initial_point=initial_point)

        M = magnetization(nqubits)
        H = ham_generator(nqubits, h, J)
        result = vqe.compute_minimum_eigenvalue(H, aux_operators=[M])

        print("Experiment complete.".ljust(30))
        print(f"GS: {result.optimal_value}")
        print(f"Magn:{result.aux_operators_evaluated}")
        session.close()

    print('noisy result for the gs is: ', result.eigenvalue)

    tmp_title = "flag" + "_" + str(J) + "_" + str(h) + "_" + "_" + str(nqubits) + ".txt"
    values = open(tmp_title, "w+")
    print(result.optimal_value, file=values)
    print((result.aux_operators_evaluated)[0][0], file=values)

    values.close()

    conv_title = str(flag)+"_"+str(J)+"_"+str(h) + "_" + "_" + str(nqubits) + "_"  + str(nreps) + "_"+str(niter) + "_"+"conv" + ".txt"
	
    conv = open(conv_title,"w+")
    print(log.values, file = conv)
	
    conv.close()
    
elif flag == 2:
    print('real hardware')

    # REAL HW
    # real device
    with Session(service=service, backend=backend) as session:

        # Set options to include the noise model
        options = Options()

        # Set number of shots, optimization_level and resilience_level
        options.execution.shots = 1024
        options.optimization_level = 3
        options.resilience_level = 1

        log = VQELog([], [])
        vqe = VQE(Estimator(session=session, options=options),
                  ansatz, optimizer, callback=log.update, initial_point=initial_point)
        M = magnetization(nqubits)
        H = ham_generator(nqubits, h, J)
        result = vqe.compute_minimum_eigenvalue(H, aux_operators=[M])

        # print("Experiment complete.".ljust(30))
        # print(f"GS: {result.optimal_value}")
        # print(f"Magn:{result.aux_operators_evaluated}")

        tmp_title = str(flag) + "_" + str(J) + "_" + str(h) + "_" + "_" + str(nqubits) + ".txt"
        values = open(tmp_title, "w+")
        print(result.optimal_value, file=values)
        print((result.aux_operators_evaluated)[0][0], file=values)

        values.close()
        conv_title = str(flag)+"_"+str(J)+"_"+str(h) + "_" + "_" + str(nqubits) + "_"  + str(nreps) + "_"+str(niter) + "_"+"conv" + ".txt"
	
        conv = open(conv_title,"w+")
        print(log.values, file = conv)
	
        conv.close()



