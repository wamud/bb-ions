# This file generates BB code circuits with the syndrome extraction circuits as usual (i.e. following Edwin Tham, Min Ye ... Delfosse "modules" paper's Algorithm 2) but with the noise model based on Helios' quantum computer.


import sys
import os
import deltakit_stim

# Save deltakit_stim as stim so any other calls to import stim (including from other packages) just use deltakit_stim:
injections = ['stim._detect_machine_architecture', 'stim._stim_polyfill', 'stim']
for namespace in injections:
    sys.modules[namespace] = sys.modules["deltakit_{namespace}".format(namespace=namespace)]


sys.path.append(os.path.abspath("../src"))
from bb_ions import *


memory_basis = 'X'
noise = 'helios'
num_syndrome_extraction_cycles = 10
seq_gates = True
exclude_opp_basis_detectors = True
leakage = True
only_CZs = True
leakage_repumping = True
cycles = 4
# p = 0.001

for p in [0.002, 0.003, 0.004, 0.005, 0.006]:
    for code in [bb5_48_4_7(), bb5_60_4_8()]:

        filename = f"n={code.n},k={code.k},d={code.d_max},p={p},noise={noise},leakage={leakage},leakage_repumping={leakage_repumping},repumping_cycles={cycles},r={num_syndrome_extraction_cycles},seq_gates={seq_gates},b={memory_basis},excl_opp_b_detectors={exclude_opp_basis_detectors},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
        
        print(filename)

        circuit = make_BB_circuit(  # see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs
            code,
            p,  
            errors = helios_errors(p),
            idle_during = helios_idle_errors(p),
            num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,  
            memory_basis = memory_basis,
            sequential_gates = seq_gates, 
            exclude_opposite_basis_detectors = exclude_opp_basis_detectors,
            reuse_check_qubits = True,
            leakage = leakage,
            only_CZs = only_CZs,
            leakage_repumping = leakage_repumping,
            num_repumping_cycles = cycles,
        )

        # Save circuit:
        circuit.to_file(f"../circuits/with_leakage/helios/actual_helios/higher_p/{filename}.stim")