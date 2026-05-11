# This file generates BB code circuits with the syndrome extraction circuits as usual (i.e. following Edwin Tham, Min Ye ... Delfosse "modules" paper's Algorithm 2) but with the noise model based on Helios' quantum computer.



import sys
import os
import deltakit_stim

# Save deltakit_stim as stim so any other calls to import stim (including from other packages) just use deltakit_stim:
injections = ['stim._detect_machine_architecture',
'stim._stim_polyfill',
'stim']
for namespace in injections:
    sys.modules[namespace] = sys.modules["deltakit_{namespace}".format(namespace=namespace)]

sys.path.append(os.path.abspath("../src"))
from bb_ions import *


ps = [0.001, 0.002]
seq_gates = False 
exclude_opp_basis_detectors = True  # If set to false then it includes detectors on X (Z) stabiliser measurement results during Memory Z (X) -- i.e. allows correlated decoding
swap_LRC = False


for code in [bb6_248_10_14()]:
    print(f"[[{code.n}, {code.k}, {code.d_max}]]")
    
    for memory_basis in 'X': # Helios suffers dephasing idling and the CZ gates are dominated by IZ and ZI errors so do mem X to see worst-case.

        num_syndrome_extraction_cycles = code.d_max 
        noise = 'helios'  # This is purely for the filename -- make sure the 'errors' and 'idle_during' correspond

        for p in ps:

            circuit = make_BB_circuit(  # (see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs)
                code,  
                p,  
                memory_basis = memory_basis, 
                num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,
                errors = helios_errors(p),
                idle_during = helios_idle_errors(p),
                sequential_gates = seq_gates, 
                exclude_opposite_basis_detectors = exclude_opp_basis_detectors,
                only_CZs = True,
                reuse_check_qubits = True,
                swap_LRC = swap_LRC
            )

            # Save circuit:

            filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},seq_gates={seq_gates},b={memory_basis},excl_opp_b_detectors={exclude_opp_basis_detectors},swap_LRC={swap_LRC},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
            
            circuit.to_file(f"../circuits/helios/parallel_gates/{filename}.stim")



