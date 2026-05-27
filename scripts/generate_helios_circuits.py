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


code = gross_code()
memory_basis = 'X'

noise = 'helios'
loss = False
num_syndrome_extraction_cycles = 12
seq_gates = True
exclude_opp_basis_detectors = True
leakage_heralds = True

loss_heralds = False

for p in [0.001, 0.002]:
    for swap_LRC in [True, False]:
        for leakage in [True, False]:

            filename = f"n={code.n},k={code.k},d={code.d_max},p={p},noise={noise},leakage={leakage},leakage_heralds={leakage_heralds},loss={loss},loss_heralds={loss_heralds},r={num_syndrome_extraction_cycles},seq_gates={seq_gates},b={memory_basis},excl_opp_b_detectors={exclude_opp_basis_detectors},swap_LRC={swap_LRC},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
            
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
                leakage_heralds = leakage_heralds,
                swap_LRC = swap_LRC,
                only_CZs = True
            )

            # svg = circuit.without_noise().diagram('timeline-svg') # diagram without noise
            # display(svg)
            # svg_string = str(svg)
            # with open(f"scrap.svg", "w", encoding="utf-8") as f: f.write(svg_string)

            # Save circuit:
            circuit.to_file(f"../circuits/with_leakage/helios/{filename}.stim")