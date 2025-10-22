import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *

# Manually define a code:
# # E.g. [[30, 4, 5]] BB5 (weight-5 stabilisers) code from Ye Delfosse long chain [2503.22071], Table II
# l = 5
# m = 3
# # A = x^0 + x
# # B = x^0 + y + x^2*y^2
# Aij = [(0, 0), (1, 0)]          # the powers (i, j) of each term x^i * y^j in A
# Bij = [(0, 0), (0, 1), (2, 2)]
# code = get_code_params(l, m, Aij, Bij)

# Use a predefined function from bbparamfuncs:
code = gross_code() 
# code = bb6_108_code()
# code = bb5_120_8_8_code()
# code = bb6_90_8_10_code()



# Options:
memory_basis = 'Z' 
num_syndrome_extraction_cycles = code.d_max # original BB paper used d
ps = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006]

noise = 'tham_modules' 


# Generate circuits:
for p in ps:
    
    if 'tham_modules' in noise:
        errors = tham_modules_errors(p)
        idle_during = tham_modules_idle_errors(p)
    # if 'our_modules' in noise:

    if 'zero_idling' in noise:
        idle_during = zero_idling()

    circuit = make_circuit(  # (see src/bb_ions/circfuncs for explanation of make_circuit inputs)
        code,  
        memory_basis,  
        p,  
        num_syndrome_extraction_cycles,  
        errors,
        idle_during,
        sequential_gates = False, 
        exclude_opposite_basis_detectors = True,
        reuse_check_qubits = True,  
    )

    # Save circuit:
    filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},b={memory_basis},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
    circuit.to_file(f"../circuits/{filename}.stim")
    
    #  svg = str(circuit.diagram("timeline-svg"))
    # with open(f"example_circuit_diagrams/{filename}.svg", "w", encoding="utf-8") as f: f.write(svg)
