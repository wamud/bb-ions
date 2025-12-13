import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *





ps = [0.002, 0.003, 0.004, 0.005, 0.006]
seq_gates = True
exclude_opposite_basis_detectors = True  # If set to false then it includes detectors on X (Z) stabiliser measurement results during Memory Z (X) -- i.e. allows correlated decoding


noise = 'helios'

# Generate circuits:
for code in [bb6_72_12_6_code(), bb6_90_8_10_code(), bb6_108_code(), gross_code(), two_gross_code()]:
    
    num_syndrome_extraction_cycles = code.d_max
    memory_basis = 'X'  # Helios suffers dephasing idling and the CZ gates are dominated by IZ and ZI errors so do mem X to see worst-case.

    for p in ps:

        circuit = make_BB_circuit(  # (see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs)
            code,  
            p,  
            memory_basis = memory_basis, 
            num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,
            errors = helios_errors(p),
            idle_during = helios_idle_errors(),
            sequential_gates = seq_gates, 
            exclude_opposite_basis_detectors = exclude_opposite_basis_detectors,
            only_CZs = True, # NECESSARY FOR HELIOS 
        )

        # Save circuit:

        prefix = "include" if exclude_opposite_basis_detectors == False else "exclude"
        
        #filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},T2={T2},pause={pause_time},seq_gates={seq_gates},b={memory_basis},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
        
        filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},seq_gates={seq_gates},b={memory_basis},excl_opp_b_detectors={exclude_opposite_basis_detectors},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
        
        circuit.to_file(f"../circuits/{noise}/{prefix}_opp_basis_detectors/{filename}.stim")

