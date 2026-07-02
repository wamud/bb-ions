# This file generates BB code circuits with the syndrome extraction circuits as usual (i.e. following Edwin Tham, Min Ye ... Delfosse "modules" paper's Algorithm 2) but with the noise model based on Helios' quantum computer.


import sys
import os
import deltakit_stim

sys.modules["stim"] = deltakit_stim 

sys.path.append(os.path.abspath("../src"))
from bb_ions import *


# code = bb5_60_4_8()

noise = 'helios'
memory_basis = 'X'
num_syndrome_extraction_cycles = 18
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
p = 0.001



for code in [bb6_756_code()]:
    filename = f"n={code.n},k={code.k},d={code.d_max},p={p},noise={noise},leakage={leakage},loss={loss},leakage_heralds={leakage_heralds},leakage_repumping={leakage_repumping},repumping_cycles={cycles},swapLRC={swapLRC},rounds={num_syndrome_extraction_cycles},seq_ops={seq_ops},seq_meas={seq_meas},b={memory_basis},excl_opp_detectors={exclude_opp_basis_detectors},l={code.l},m={code.m},A={'+'.join(f'x{i}y{j}' for i, j in code.Aij)},B={'+'.join(f'x{i}y{j}' for i, j in code.Bij)}"
    
    print("Creating: ",filename)

    circuit = make_BB_circuit(  # see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs
        code,
        errors = helios_errors(p),
        idle_during = helios_idle_errors(p),
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
        leakage_heralds=leakage_heralds
    )

    # Save circuit:
    circuit.to_file(f"../circuits/leakage_and_loss/all_codes/{filename}.stim")
