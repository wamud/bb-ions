import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *





ps = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
seq_gates = True
noise = 'uniform_plus_shift_and_shuttle_w_dephasing_idling'

# Generate circuits:
for code in [bb6_90_8_10_code(), bb6_108_code(), gross_code(), two_gross_code()]:
    num_syndrome_extraction_cycles = code.d_max

    for p in ps:
        for memory_basis in 'XZ':
            for T2 in [1, 5, 10, 50]:
                for pause_time in [0]:

                    errors = our_uniform_plus_shift_and_shuttle(p, T2)
                    idle_during = dephasing_idle_errors(T2)

                    idle_during["pause"].p = p_idle_dephasing(pause_time, T2)

                    circuit = make_BB_circuit(  # (see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs)
                        code,  
                        p,  
                        memory_basis = memory_basis,
                        num_syndrome_extraction_cycles = num_syndrome_extraction_cycles,
                        errors = errors,
                        idle_during = idle_during,
                        sequential_gates = seq_gates, 
                    )

                    # Save circuit:
                    filename = f"nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},noise={noise},r={num_syndrome_extraction_cycles},T2={T2},pause={pause_time},seq_gates={seq_gates},b={memory_basis},l={code.l},m={code.m},A='{''.join(str(x) + str(y) for x, y in code.Aij)}',B='{''.join(str(x) + str(y) for x, y in code.Bij)}'"
                    circuit.to_file(f"../circuits/{noise}/{filename}.stim")
                    
                    # svg_string = str(circuit.diagram("timeline-svg"))
                    # with open(f"scrap.svg", "w", encoding="utf-8") as f: f.write(svg_string)

