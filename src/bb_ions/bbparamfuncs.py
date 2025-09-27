'''Once the parity check matrices of a BB code have been constructed using the functions in bbfuncs.py, it is useful to find the paramaters of the code.
I.e. how many logical qubits it contains, what each logical qubit's logical operators (X_L and Z_L) are and the distance of the code.

Note on code distance:
Code distance is the smallest number of single-qubit operations required to cause a logical error. In the case of this being different for different
logical qubits in the BB code, the smallest distance is given. Note however that the distance returned from the bposd simulations is a 'maximum' distance.
That is, while the smallest layout of single-qubit operations that it found that caused a logical error might be, for example, d = 5, it can not guarantee
that a layout with a smaller number, for example d = 4, doesn't exist. Hence d = 5 is the maximum possible distance that the code could have. The more 
trials that are run the more confident we become that this is its actual distance however it's not certain.
'''

import numpy as np
import bposd
from bposd import css
from autqec.utils.qec import *


'''get_code_params
Given Hx and Hz, the X and Z parity check matrices of a BB code, this function used Joschka Roff's bposd package to return the number of data qubits n and number of logical qubits k
'''
def get_code_params(Hx, Hz):
    code = css.css_code(hx = Hx, hz = Hz)
    return n, k


'''bposd_logical_ops --(Recommended to use autqec_logical_ops function instead to return a canonical set of logical ops)
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

'''autqec_logical_ops
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








