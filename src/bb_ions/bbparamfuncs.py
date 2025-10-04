''' bbparamfuncs
Once the parity check matrices of a BB code have been constructed using the functions in bbfuncs.py, it is useful to find the paramaters of the code.
I.e. how many logical qubits it contains, what each logical qubit's logical operators (X_L and Z_L) are and the distance of the code.

Note on code distance:
Code distance is the smallest number of single-qubit operations required to cause a logical error. In the case of this being different for different
logical qubits in the BB code, the smallest distance is given. Note however that the distance returned from the bposd simulations is a 'maximum' distance.
That is, while the smallest layout of single-qubit operations that it found that caused a logical error might be, for example, d = 5, it can not guarantee
that a layout with a smaller number, for example d = 4, doesn't exist. Hence d = 5 is the maximum possible distance that the code could have. The more 
trials that are run the more confident we become that this is its actual distance however it's not certain.
'''

import numpy as np
import json
import bposd
from bposd import css
from bposd import css_decode_sim
from autqec.utils.qec import *
from .bbfuncs import *
from .modulefuncs import *

mpow = np.linalg.matrix_power

class Code:
    def __init__(self, l, m, Aij, Bij, ATij, BTij, A, B, Hx, Hz, Lx, Lz, d_max, n, k, Junion, JTunion):
        self.l = l
        self.m = m
        self.Aij = Aij
        self.Bij = Bij
        self.ATij = ATij
        self.BTij = BTij
        self.A = A
        self.B = B
        self.Hx = Hx
        self.Hz = Hz
        self.Lx = Lx
        self.Lz = Lz
        self.d_max = d_max
        self.n = n
        self.k = k
        self.Junion = Junion
        self.JTunion = JTunion




''' find_d_max
Given the X and Z parit check matrices of a css code, this function finds the maximum distance it could have. I.e. it runs simulations with bposd decoder to find the distance of the code. We call this maximum distance as there could be a logical operator with a smaller weight that would cause a logical error (i.e., a smaller distance) but the simulations didn't see it. Increase target_runs for more certainty in d_max'''
def find_d_max(Hx, Hz):
    
    osd_options={
    'target_runs': 2000, 'xyz_error_bias': [1, 1, 1], 'bp_method': "minimum_sum", 'ms_scaling_factor': 0.05, 'osd_method': "osd_cs", 'osd_order': 4, 'channel_update': None, 'seed': 42, 'max_iter': 9, # 'output_file': "test.json", 
    'error_bar_precision_cutoff': 1e-6
    }

    error_rate = 0.1
    bb5 = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = error_rate, **osd_options) 

    results = json.loads(bb5.output_dict())
    d_max = results["min_logical_weight"]

    return d_max



''' get_code_params
Using the l, m, A indices (Aij) and B incdices (Bij) of a bicycle bivariate code this function constructs a code class which contains its parity check matrices, logical operators, d_max (could be a smaller distance but the simulations found d_max) n, k l, m, A and B.
Note that Aij and Bij are inserted are tuples (i, j) for each of the terms x^i y^j in A and B. 
E.g. if
A = x^0 + x
B = x^0 + y + x^2*y^2   ([[30, 4, 5]] code from Table II of [2503.22071])
Then
Aij = [(0, 0), (1, 0)]
Bij = [(0, 0), (0, 1), (2, 2)]
The A and B matrices returned by this function are the actual matrices A, B
'''
def get_code_params(l, m, Aij, Bij):

    # Sorting indices into I(A), I(B), J(A), J(B):
    IA, JA = findIJ(Aij)
    IB, JB = findIJ(Bij)
    Junion = sorted(set(JA + JB))

    # Sorting indices into I(A^T), I(B^T), J(A^T), J(B^T):
    ATij = [(-i, -j) for i, j in Aij]
    BTij = [(-i, -j) for i, j in Bij]
    IAT, JAT = findIJ(ATij)
    IBT, JBT = findIJ(BTij)
    JTunion = sorted(set(JAT + JBT))


    # Num qubits:
    n = 2 * l * m

    # Define matrices A and B:    
    ## Setup:
    A = 0; B = 0
    x = make_x(l, m)
    y = make_y(l, m)
    
    ## A = A1 + A2 + ...
    for (i, j) in Aij: 
        A = mpow(x, i) @ mpow(y, j) + A
    
    ## B = B1 + B2 + ...
    for (i, j) in Bij:
        B = mpow(x, i) @ mpow(y, j) + B

    # Parity check matrices:
    Hx = np.hstack((A, B))
    Hz = np.hstack((B.T, A.T))

    # Logical operators:
    Lx, Lz = autqec_logical_ops(Hx, Hz)

    # Num. logical qubits:
    k = len(Lx)

    # Maximum distance:
    d_max = 5      #find_d_max(Hx, Hz)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -- hard-coded while working on 30,4,5 

    code = Code(l, m, Aij, Bij, ATij, BTij, A, B, Hx, Hz, Lx, Lz, d_max, n, k, Junion, JTunion)

    return code




''' bposd_logical_ops --(Recommended to use autqec_logical_ops function instead to return a canonical set of logical ops)
Given the parity check matrices of a BB code, let's get the logical operators for each of its logical qubits.
This function does this using Joschka Roffe's bposd package.
Checking the anticommutation relations of the logical operators is done with anticommute_matrix = Lx @ Lz^T % 2. It checks each Lx against every Lz. The anticommute matrix shows that bposd returns a generating set of logical operators that is not canonical. I.e. the anticommute matrix has rank = k but is not I_k. Row-reducing the anticommute matrix until it is I_k (and performing the corresponding multiplication of logical operators) will give a canonical set where each logical operator anticommutes with its partner on the same logical qubit but commutes with all the others.
This is done within the autqec_logical_ops function so it may be used instead.
'''
def bposd_logical_operators(Hx, Hz):
    code = css.css_code(hx = Hx, hz = Hz)
    # n = code.N
    # k = code.K
    # d = code.D
    # print(f'[[{n}, {k}, {d}]]\n')

    # Look at logical ops:
    Lx = code.lx.toarray()
    Lz = code.lz.toarray()

    # Check anticommutations between logical operators (a 1 implies anticommutes)
    anticommute_matrix = Lx @ Lz.T % 2
    rank = np.linalg.matrix_rank(anticommute_matrix)
    # Note if each logical operator anticommutes with exactly one other then you get the identity matrix out
    # You might also get an equivalent binary matrix out, i.e. of the same rank. If its rank is equal to number of logical qubits then you can multiply logical operators together to eventually just get pairs that anticommute, i.e. turn the matrix into the identity matrix)
    assert(rank == code.K)
        # Correct anticommutation relations between logical operators

    return Lx, Lz

''' autqec_logical_ops
Given the parity check matrices of a BB code, let's get the logical operators for each of its logical qubits.
This function does this using Hasan Sayginel's autqec package.
Checking the anticommutation relations of the logical operators is done with anticommute_matrix = Lx @ Lz^T % 2.
It checks each Lx against every Lz.
We check that the anticommmute matrix is exactly the identity, implying we have a canonical set of logial operators, where each logical operator anticommutes with its partner on the same logical qubit but commutes with all the others.
Note: if n and k have already been found they can be fed to this function. Alternatively it will find n and k using Joschka Roffe's bposd package.
'''
def autqec_logical_ops(Hx, Hz, n = None, k = None):
    
    if n is None or k is None:
        code = css.css_code(hx = Hx, hz = Hz)
        n = code.N
        k = code.K
    
    ''' 
    First add zero matrices to convert parity check matrices to symplectic form for autqec. 
    e.g. XYZIZ in symplectic form 
    = [X part | Z part]
    = [X X I I I | I Z Z I Z]
    = [1 1 0 0 0 | 0 1 1 0 1]

    Now because our BB code is CSS the combined parity check matrix is simply
     [  Hx  |  0 
        0   | Hz  ]
    '''

    zeros = np.zeros_like(Hx)
    H_symp = np.array(np.vstack((np.hstack((Hx,zeros)),np.hstack((zeros,Hz)))),dtype=int)

    # Row reduce and find logical operators:
    H_symp_rref, _, transform_rows, transform_cols = rref_mod2(H_symp)
    H_symp_rref = H_symp_rref[~np.all(H_symp_rref == 0, axis=1)]
    H_symp_rref_og_basis = H_symp_rref@inv_mod2(transform_cols)
    assert H_symp_rref_og_basis.shape[0] == n - k
    assert H_symp_rref_og_basis.shape[1] == 2 * n

    G, LX_symplectic, LZ_symplectic, D = compute_standard_form(H_symp_rref_og_basis)

    # We have a CSS code so cut off the empty Z and X parts of the symplectic representations of Lx and Lz:
    Lx = LX_symplectic[:, :n]
    Lz = LZ_symplectic[:, n:]

    # Verify anticommutation of logical operators and that they're canonical:
    anticommute_matrix = (Lx @ Lz.T) % 2
    ident_matrix = np.eye(Lx.shape[0])
    assert(np.array_equal(anticommute_matrix, ident_matrix))

    return Lx, Lz








