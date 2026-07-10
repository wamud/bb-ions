import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF


# # Uncomment code(s) of choice:
codes = [
    bb5_30_4_5_code(), 
    bb5_48_4_7_code(), 
    bb6_72_12_6_code(), 
    bb6_90_8_10_code(), 
    bb6_108_code(), 
    gross_code(), 
    two_gross_code() 
]


p = 0.01
memory_basis = 'Z'

code = codes[int(sys.argv[1])]


msg = f"- [[{code.n}, {code.k}, {code.d_max}]]"
print(msg)


circuit = make_BB_circuit(  # see src/bb_ions/circfuncs for explanation of make_BB_circuit inputs
    code,  
    p,  
    errors = tham_modules_errors(p),
    idle_during = tham_modules_idle_errors(p),
    num_syndrome_extraction_cycles = 4,  
    memory_basis = memory_basis,
    sequential_gates = True, 
    exclude_opposite_basis_detectors = True,
    reuse_check_qubits = True,  
)

# # Generate maxSAT problem of the circuit's distance, that other tools can solve.
# with open("problem.wcnf", "w") as file:
    # file.write(circuit.shortest_error_sat_problem())

# # Solve the problem using python-sat (slower but compatible with Apple / ARM )
# wcnf = WCNF(from_file="problem.wcnf")


wcnf = WCNF(from_string = circuit.shortest_error_sat_problem())

with RC2(wcnf) as rc2:
    msg1 = rc2.compute()
    msg2 = str(rc2.cost)

    print(f"     {msg} : {msg2}")

    with open("circ_distance_distance.txt", "a") as file:
        file.write(f"{msg} : {msg2}\n")
    with open("circ_distance_solution.txt", "a") as file:
        file.write(f"{msg}\n   {msg1}\n")

