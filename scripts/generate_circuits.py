import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *



# Use a predefined function from bbparamfuncs to define code:
# code = gross_code() 
# code = bb6_108_code()
# code = bb5_120_8_8_code()
# code = bb6_90_8_10_code()



ps = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006]

noise = 'tham_modules' 


# Generate circuits:
for code in [bb6_90_8_10_code(), bb6_108_code, gross_code()]:
    for p in ps:

        num_syndrome_extraction_cycles = code.d_max
        sequential_gates = True
        memory_basis = 'Z'

        errors = tham_modules_errors(p)
        idle_during = tham_modules_idle_errors(p)

        circuit = make_BB_circuit(  # (see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs)
            code,  
            p,  
            memory_basis = 'Z',
            num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,
            errors = errors,
            idle_during = idle_during,
            sequential_gates = sequential_gates, 
            exclude_opposite_basis_detectors = True,
            reuse_check_qubits = True,  
        )

        # Save circuit:
        filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},seq_gates={sequential_gates},b={memory_basis},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
        circuit.to_file(f"../circuits/{noise}_noise/{filename}.stim")
        

# svg_string = str(circuit.diagram("timeline-svg"))
# with open(f"scrap.svg", "w", encoding="utf-8") as f: f.write(svg_string)