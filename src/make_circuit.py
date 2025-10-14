from bb_ions import *
import stim


# [[30, 4, 5]] From Ye Delfosse long chain [2503.22071] Table II
l = 5
m = 3
# A = x^0 + x
# B = x^0 + y + x^2*y^2
Aij = [(0, 0), (1, 0)]          # the powers (i, j) of x^i * y^j
Bij = [(0, 0), (0, 1), (2, 2)]
d = 5

# [[48, 4, 7]] from Ye Delfosse long chain [2503.22071] Table II
#l = 8
#m = 3
## A = x^0 + x
## B = x^0 + y + x^3 * y^2
#Aij = [(0, 0), (1, 0)]
#Bij = [(0, 0), (0, 1), (3, 2)]
#d = 7

# # [[108, 8, 10]] from BB paper [2308.07915] Table III
# l = 9
# m = 6
# # A = x^3 + y + y^2
# # B = y^3 + x + x^2
# Aij = [(3, 0), (0, 1), (0, 2)]
# Bij = [(0, 3), (1, 0), (2, 0)]
# d = 10


# # [[144, 12, 12]] 'gross code' from BB paper [2308.07915] Table III
# l = 12
# m = 6
# # A = x^3 + y + y^2
# # B = y^3 + x + x^2
# Aij = [(3, 0), (0, 1), (0, 2)]
# Bij = [(0, 3), (1, 0), (2, 0)]
# d = 12


# # [[288, 12, 18]] 'two gross'  from BB paper [2308.07915] Table III
# l = 12
# m = 12
# # A = x^3 + y^2 + y^7
# # B = y^3 + x^1 + x^2
# Aij = [(3, 0), (0, 2), (0, 7)]
# Bij = [(0, 3), (1, 0), (2, 0)]
# d = 18

code = get_code_params(l, m, Aij, Bij, d)
memory_basis = 'Z' # preserve logical 0 or + if basis is Z or X
p = 0.001
noise = 'longchain' # current options are 'longchain' or 'uniform'
num_syndrome_extraction_cycles = code.d_max


circuit = make_circuit(
    code,
    memory_basis,
    p,
    noise,
    num_syndrome_extraction_cycles,
    sequential = True,
    exclude_opposite_basis_detectors = False,
    reuse_check_qubits = True
)

# circuit.to_file(f"../circuits/nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},b={memory_basis},noise={noise},r={num_syndrome_extraction_cycles},code=BB,l={l},m={m},A='{''.join(str(x) + str(y) for x, y in Aij)}',B='{''.join(str(x) + str(y) for x, y in Bij)}'.stim")


circuit.to_file(f"../circuits/hi.stim")
