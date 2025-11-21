import stim
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF


gross_code = stim.Circuit.from_file("../circuits/relaybp_circuits/circuit=bicycle_bivariate_144_12_12_memory_Z,distance=12,rounds=12,error_rate=0.001,noise_model=uniform_circuit,basis=CX,A=x^3+y+y^2,B=y^3+x+x^2.stim")

code_72 = stim.Circuit.from_file("../circuits/relaybp_circuits/circuit=bicycle_bivariate_72_12_6_memory_Z,distance=6,rounds=6,error_rate=0.001,noise_model=uniform_circuit,basis=CX,A=x^3+y+y^2,B=y^3+x+x^2.stim")

circuits = [code_72, gross_code]

circuit = circuits[int(sys.argv[1])]

if int(sys.argv[1]) == 0:
    n = 72
    k = 12 
    d_max = 6
elif int(sys.argv[1]) == 1:
    n = 144
    k = 12
    d_max = 12
else:
    raise ValueError("usage: python3 find_relay_circuit_level_distance.py 0/1")

msg = f"- [[{n}, {k}, {d_max}]]"
print(msg)

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

    with open("relay_circ_distance_distance.txt", "a") as file:
        file.write(f"{msg} : {msg2}\n")
    with open("relay_circ_distance_solution.txt", "a") as file:
        file.write(f"{msg}\n   {msg1}\n")

