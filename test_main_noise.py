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


def old_ham_generator(N, h, J):
    H = 0

    # field term
    for i in range(N):
        H += h * ((I ^ (i)) ^ (Z) ^ (I ^ (N - i - 1)))
    # interaction term
    for k in range(N - 1):
        H += J * ((I ^ (k)) ^ (X ^ X) ^ (I ^ (N - k - 2)))
    # PBC
    H += J * ((X) ^ (I ^ (N - 2)) ^ (X))
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


# definition magnetization
def old_magnetization(N):
    M = 0

    # field term
    for i in range(N):
        M += ((I ^ (i)) ^ (X) ^ (I ^ (N - i - 1)))
    return M / N

nqubits=5 #6
h = 0.2 #or -0.2 or 0
J=1
H = ham_generator(nqubits, h, J)
M = magnetization(nqubits)
#H, M


sol = NumPyMinimumEigensolver().compute_minimum_eigenvalue(H, aux_operators=[M])

print('correct solution is: ', sol.eigenvalue)

service = QiskitRuntimeService(channel='ibm_quantum')
sim = "ibmq_qasm_simulator"

ansatz = RealAmplitudes(num_qubits=nqubits, reps=3)
#ansatz = EfficientSU2(nqubits, reps=3, entanglement='linear', insert_barriers=False)
#ansatz_true = EfficientSU2(nqubits, reps=3, entanglement='linear', insert_barriers=True)

optimizer = SPSA(maxiter=200)

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

# Make a noise model

fake_backend = FakeAuckland()
'''noise_model = NoiseModel.from_backend(fake_backend)
basis_gates=fake_backend.configuration().basis_gates
coupling_map=fake_backend.configuration().coupling_map
'''

noise_model = NoiseModel.from_backend(backend)
basis_gates=backend.configuration().basis_gates
coupling_map=backend.configuration().coupling_map

# Set options to include the noise model
options = Options()
options.simulator = {
    "noise_model": noise_model,
    "basis_gates":basis_gates,
    "coupling_map": coupling_map,
    "seed_simulator": 42
}

# Set number of shots, optimization_level and resilience_level
options.execution.shots = 1000
options.optimization_level = 3
options.resilience_level = 0

#noisy case
with Session(service=service, backend=sim) as session:
    estimator=Estimator(session=session, options=options)
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




