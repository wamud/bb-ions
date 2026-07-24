# This file generates BB code circuits with the syndrome extraction circuits as usual (i.e. following Edwin Tham, Min Ye ... Delfosse "modules" paper's Algorithm 2) but with the noise model based on Helios' quantum computer.
import sys
import os
import deltakit_stim
sys.modules["stim"] = deltakit_stim 
sys.path.append(os.path.abspath("../src"))
from bb_ions import *
import json

memory_basis = 'X' # mem X is worst-case in Helios due to two-qubit gate and idling being dominate by Z errors
num_syndrome_extraction_cycles = 20
only_CZs = True
seq_ops = True
seq_meas = False
leakage = False
loss = False
leakage_repumping = False
cycles = 3
swapLRC = False
leakage_heralds = False





codes = [interchanged_bb5_60_4_8(), bb6_72_12_6(), gross_code()]

for exclude in [True, False]:
    for code in codes:
        for p in [0.001, 0.002, 0.003]:


            # if seq_ops == True and seq_meas == True:
            #     max_parallel_1q_ops = 16 * code.m
            #     max_parallel_2q_ops = 4 * code.m
            # else :
            #     max_parallel_1q_ops = np.inf
            #     max_parallel_2q_ops = np.inf


            max_parallel_1q_ops = np.inf
            max_parallel_2q_ops = code.m




            filename = f"n={code.n},k={code.k},d={code.d_max},l={code.l},m={code.m},A={'+'.join(f'x{i}y{j}' for i, j in code.Aij)},B={'+'.join(f'x{i}y{j}' for i, j in code.Bij)},p={p},leakage={leakage},loss={loss},leak_heralds={leakage_heralds},leak_repump={leakage_repumping},repump_cycles={cycles},swapLRC={swapLRC},r={num_syndrome_extraction_cycles},seq_ops={seq_ops},seq_meas={seq_meas},b={memory_basis},excl_opp_detectors={exclude},prllel_1q={max_parallel_1q_ops},prllel_2q={max_parallel_2q_ops}"
            
            
            print("Creating: ",filename)

            circuit = make_BB_circuit(  # see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs
                code,
                errors = helios_errors(p, code),
                idle_during = helios_idle_errors(p, code),
                num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,  
                memory_basis = memory_basis,
                sequential_operations = seq_ops, 
                sequential_measurements = seq_meas,
                exclude_opposite_basis_detectors = exclude,
                reuse_check_qubits = True,
                leakage = leakage,
                only_CZs = only_CZs,
                leakage_repumping = leakage_repumping,
                num_repumping_cycles=cycles,
                swapLRC = swapLRC,
                loss = loss,
                leakage_heralds=leakage_heralds,
                num_parallel_1q_ops = max_parallel_1q_ops, 
                num_parallel_2q_ops = max_parallel_2q_ops,
            )


            # Save circuit:
            circuit.to_file(f"../circuits/without_leakage/exclude_vs_include/{filename}.stim")
            # svg = circuit.diagram('timeline-svg')
            # svg_string = str(svg)
            # with open(f"../scrap.svg", "w", encoding="utf-8") as f: f.write(svg_string)
