# importing libraries 
import sys
import numpy as np
import matplotlib.pyplot as plt

# standard Qiskit libraries
from qiskit import QuantumCircuit, transpile, Aer, IBMQ
from qiskit.tools.jupyter import *
from qiskit.visualization import *
from ibm_quantum_widgets import *
from qiskit.providers.aer import QasmSimulator
from qiskit.opflow.primitive_ops import PauliSumOp
from qiskit.providers.aer import AerSimulator
from qiskit.opflow import Z,I,X
from qiskit.circuit.library import EfficientSU2
from qiskit.algorithms.optimizers import SPSA,  COBYLA
from qiskit.providers.basicaer import QasmSimulatorPy  # local simulator 
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.algorithms.minimum_eigensolvers import VQE
import matplotlib.pyplot as plt
from qiskit.providers.fake_provider import FakeGeneva, FakeVigo, FakeNairobi, FakeAuckland
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime import QiskitRuntimeService, Estimator, Session, Options

# useful functions
def ham_generator(N,h,J):
    H=0
    
    # field term
    for i in range (N):
        H+= h*((I^(i))^(Z)^(I^(N-i-1))) 
    # interaction term
    for k in range(N-1):
        H+= J*((I^(k))^(X^X)^(I^(N-k-2)))
    # PBC
    H+= J*((X)^(I^(N-2))^(X))           
    return H    

#function that calculates the forward derivative
def f_derivative(a,b,en,npoints):

    step = (b-a)/(npoints-1)
    for i in range(npoints):
        if a + i*step == 0:
            print('zero found in position >>', i)
            z=i
            break
    fd = (en[z+1]-en[z])/(step)
    
    print('spacing used to calculate the derivative >> ', step)
    print('numerical derivative of E(h), h->0+ >>', fd)
    return fd


#function that calculates the backward derivative
def b_derivative(a,b,en,npoints):

    step = (b-a)/(npoints-1) 
    for i in range(npoints):
        if a + i*step == 0:
            print('zero found in position >>', i)
            z=i
            break
    bd = (en[z]-en[z-1])/(step)
    
    print('spacing used to calculate the derivative >> ', step)
    print('numerical derivative of E(h), h->0- >>', bd)
    return bd

#definition magnetization
def magnetization(N):
    M=0
    
    # field term
    for i in range(N):
        M+= ((I^(i))^(X)^(I^(N-i-1)))           
    return M/N 


flag = int(sys.argv[1])

J = float(sys.argv[2])

h = float(sys.argv[3])

nqubits = int(sys.argv[4])

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
        #if (count % 80) == 0:
            #print( parameters)
log = VQELog([],[])



service = QiskitRuntimeService()
backend = "ibmq_qasm_simulator"

ansatz = EfficientSU2(nqubits, reps=3, entanglement='linear', insert_barriers=True)

optimizer = SPSA(maxiter=2000)

initial_point = np.random.random(ansatz.num_parameters)


if flag == 0:
    
    with Session(service=service, backend=backend) as session:
    
        log = VQELog([],[])
        vqe = VQE(Estimator(),
                  ansatz=ansatz, optimizer=optimizer, initial_point=initial_point, callback=log.update)
        M=magnetization(nqubits)
        H = ham_generator(nqubits,h,J)
        result = vqe.compute_minimum_eigenvalue(H,  aux_operators=[M])
   
    #print("Experiment complete.".ljust(30))
    #print(f"GS: {result.optimal_value}")
    #print(f"Magn:{result.aux_operators_evaluated}")
        
        
        tmp_title = str(flag)+"_"+str(J)+"_"+str(h)+"_"+"_"+str(nqubits)+".txt"
        values = open(tmp_title,"w+")
        print(result.optimal_value,file=values)
        print(result.aux_operators_evaluated, file = values)
        
        values.close()

elif flag == 1:
    
    with Session(service=service, backend=backend) as session:
        # Make a noise model
        fake_backend = FakeAuckland()
        noise_model = NoiseModel.from_backend(fake_backend)

        # Set options to include the noise model
        options = Options()
        options.simulator = {
            "noise_model": noise_model,
            "basis_gates": fake_backend.configuration().basis_gates,
            "coupling_map": fake_backend.configuration().coupling_map,
            "seed_simulator": 42
            }

        # Set number of shots, optimization_level and resilience_level
        options.execution.shots = 2000
        options.optimization_level = 3
        options.resilience_level = 2
    
        log = VQELog([],[])
        vqe = VQE( Estimator(session=session,options=options),
                  ansatz, optimizer, callback=log.update, initial_point=initial_point)
        M=magnetization(nqubits)
        H = ham_generator(nqubits,h,J)
        result = vqe.compute_minimum_eigenvalue(H,  aux_operators=[M])
        
   
    #print("Experiment complete.".ljust(30))
    #print(f"GS: {result.optimal_value}")
    #print(f"Magn:{result.aux_operators_evaluated}")
    
        tmp_title = str(flag)+"_"+str(J)+"_"+str(h)+"_"+"_"+str(nqubits)+".txt"
        values = open(tmp_title,"w+")
        print(result.optimal_value,file=values)
        print(result.aux_operators_evaluated, file = values)

        values.close()
