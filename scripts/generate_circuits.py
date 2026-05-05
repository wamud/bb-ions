import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *





ps = [0.001]

seq_gates = False
exclude_opposite_basis_detectors = False  # If set to false then it includes detectors on X (Z) stabiliser measurement results during Memory Z (X) -- i.e. allows correlated decoding
swapLRC = False

# Generate circuits:
for code in [two_gross_code()]:
    
    num_syndrome_extraction_cycles = code.d_max

    for p in ps:
        for memory_basis in 'X':

            noise = 'uniform'
            errors = uniform_errors(p)
            idle_during = uniform_idling(p)

            circuit = make_BB_circuit(  # (see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs)
                code,  
                p,  
                memory_basis = memory_basis,
                num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,
                errors = errors,
                idle_during = idle_during,
                sequential_gates = seq_gates, 
                exclude_opposite_basis_detectors = exclude_opposite_basis_detectors,
                swap_LRC = swapLRC
            )


            filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},seq_gates={seq_gates},b={memory_basis},excl_opp_b_detectors={exclude_opposite_basis_detectors},swap_LRC={swapLRC},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)},SE='non-interleaved''"

            circuit.to_file(f"../circuits/{noise}/{filename}.stim")   
