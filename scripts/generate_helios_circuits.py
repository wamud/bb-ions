# This file generates BB code circuits with the syndrome extraction circuits as usual (i.e. following Edwin Tham, Min Ye ... Delfosse "modules" paper's Algorithm 2) but with the noise model based on Helios' quantum computer.
import sys
import os
import deltakit_stim
sys.modules["stim"] = deltakit_stim 
sys.path.append(os.path.abspath("../src"))
from bb_ions import *
import json

noise = 'helios'
memory_basis = 'X' # mem X is worst-case in Helios due to two-qubit gate and idling being dominate by Z errors
num_syndrome_extraction_cycles = 20
only_CZs = True
seq_ops = True
seq_meas = True
leakage = True
loss = True
leakage_repumping = True
cycles = 3
swapLRC = False
leakage_heralds = False
exclude_opp_basis_detectors = True


for p in [0.001]:
    for code in [four_twenty_code(), bb6_210_8_18()]:
        

        max_parallel_1q_ops = 16 * code.m
        max_parallel_2q_ops = 4 * code.m


        filename = f"n={code.n},k={code.k},d={code.d_max},l={code.l},m={code.m},p={p},leakage={leakage},loss={loss},leak_heralds={leakage_heralds},leak_repump={leakage_repumping},repump_cycles={cycles},swapLRC={swapLRC},r={num_syndrome_extraction_cycles},seq_ops={seq_ops},seq_meas={seq_meas},prllel_1q={max_parallel_1q_ops},prllel_2q={max_parallel_2q_ops},b={memory_basis},excl_opp_detectors={exclude_opp_basis_detectors},A={'+'.join(f'x{i}y{j}' for i, j in code.Aij)},B={'+'.join(f'x{i}y{j}' for i, j in code.Bij)}"
        
        print("Creating: ",filename)

        circuit = make_BB_circuit(  # see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs
            code,
            errors = helios_errors(p, code),
            idle_during = helios_idle_errors(p, code),
            num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,  
            memory_basis = memory_basis,
            sequential_operations = seq_ops, 
            sequential_measurements = seq_meas,
            exclude_opposite_basis_detectors = exclude_opp_basis_detectors,
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
        circuit.to_file(f"../circuits/leakage_and_loss/all_codes/{filename}.stim")
