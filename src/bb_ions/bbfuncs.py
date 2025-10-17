''' bbfuncs
A Bicycle Bivariate (BB) code [2308.07915] is a CSS code and thus has separate Hx and Hz parity check matrices. 
These functions are for constructing the parity check matrices of a BB code.'''

import numpy as np
import scipy
from ldpc import mod2
mpow = np.linalg.matrix_power

''' make_s
makes a cyclic matrix S. This is an idenity matrix with every 1 cyclically shifted to the right by one'''
def make_s(dim):
  s = np.zeros((dim,dim), dtype = int)

  for i in range(dim):
    s[i % dim, (i + 1) % dim] = 1

  return s

''' make_x
Makes the matrix x, which is x := S_l ⊗ I_m'''
def make_x(l, m):
  s_l = make_s(l)
  ident_m = np.eye(m, dtype = int)

  x = np.kron(s_l, ident_m)
  return x

''' make_y
Makes the matrix y, which is y := I_l ⊗ S_m'''
def make_y(l,m):

  ident_l = np.eye(l, dtype = int)
  s_m = make_s(m)

  y = np.kron(ident_l, s_m)
  return y

''' make_z
Makes the matrix z, which is x * y or, equivalently, S_l ⊗ S_m'''
def make_z(l,m):
  s_l = make_s(l)
  s_m = make_s(m)

  z = np.kron(s_l, s_m)

  return z

''' make_xyz
Uses above functions to make x, y and z'''
def make_xyz(l,m):
  x = make_x(l, m)
  y = make_y(l, m)
  z = make_z(l, m)

  return x, y, z

''' make_parity_check_matrices
The parity check matrices of a BB code are Hx = [A|B] and Hz = [B^T | A^T], where A and B are sums of x^i * y^j , i.e. S_l ⊗ S_m where S_a is a cyclic matrix of dimension a × a. This function takes as inputs Aij and Bij which are the powers (i, j) of each term in A and B. E.g. if A = x^0 + x^2 * y^2 then Aij = [(0, 0), (2, 2)]. The function then returns the two parity check matrices Hx and Hz.'''
def make_parity_check_matrices(l, m, Aij, Bij):
    ## Setup:
    A = 0; B = 0
    x = make_x(l, m)
    y = make_y(l, m)

    ## A = A1 + A2 + ...
    for (i, j) in Aij: 
        A = (mpow(x, i) @ mpow(y, j) + A ) % 2

    ## B = B1 + B2 + ...
    for (i, j) in Bij:
        B = ( mpow(x, i) @ mpow(y, j) + B ) % 2

    Hx = np.hstack((A, B))
    Hz = np.hstack((B.T, A.T))

    return Hx, Hz

''' make_sparse_parity_check_matrices
Same as make_parity_check_matrices (above) but returns sparse matrices. Could make this faster by working with sparse matrices from the beginning (i.e. in make_x, make_y etc.)'''
def make_sparse_parity_check_matrices(l, m, Aij, Bij):
    
    Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)
    
    Hx = scipy.sparse.csr_matrix(Hx)
    Hz = scipy.sparse.csr_matrix(Hz)

    return Hx, Hz



''' is_valid
Given two parity check matrices of a CSS code, tests they have the same number of qubits and that their stabilisers commute'''
def is_valid(Hx, Hz):

    if Hx.shape[1] != Hz.shape[1]:  # Hx and Hz must have the same number of columns
        return False

    anticom_matrix = (Hx @ Hz.T) % 2 
    if np.any(anticom_matrix): # if any 1's (indicating anticommuting stabilisers)
        return False

    return True

''' find_k
Finds the number of logical operators of a CSS code given its two parity check matrices'''
def find_k(Hx, Hz):

    n = Hx.shape[1]
    rank_Hx = mod2.rank(Hx)
    rank_Hz = mod2.rank(Hz)
    
    k = n - rank_Hx - rank_Hz
    
    return k



''' test_stabs_commute
Given X and Z parity check matrices, this function tests that all the stabilisers commute by testing Hx Hz^T = [0] mod 2
(they should have an even number of 1's in corresponding positions so give the zero matrix)'''
def test_stabs_commute(Hx, Hz):
  test = (Hx @ Hz.T) % 2 # should give all-zero matrix zeros
  if np.any(test) == True: # any 1's (indicating anticomm.)
    raise ValueError("Tragédie, tragédie ! Not all the stabilisers commute")


''' verify_ones
For each matrix in args, verify that there is just one one per row and the rest are zeros:
'''
def verify_ones(*args):
  for i, A in enumerate(args):
    sumlist = np.sum(A, axis = 1) # sum of each row
    num_zeros = np.sum(A == 0, axis = 1) # num of zeros in each row

    test1 = np.all(sumlist == 1)
    numcols = A.shape[1]
    test2 = np.all(num_zeros == numcols - 1)

    if not (test1 and test2):
      print(f"Tragédie, tragédie ! The {i}-th matrix given to verify_ones is not exactly one one per row")




# # E.g.

# # [[30, 4, 5]] From Ye Delfosse long chain [2503.22071]

# # l = 5, m = 3, A = x^0 + x , B = x^0 + y + x^2*y^2  ## (from table II)

# l = 5
# m = 3

# x, y, z = make_xyz(l,m)

# x0 = np.linalg.matrix_power(x,0)

# A1 = x0
# A2 = np.zeros((l*m, l*m), dtype = int)
# A3 = x

# B1 = x0
# B2 = y
# B3 = z @ z

# A = A1 + A2 + A3
# B = B1 + B2 + B3

# Hx = np.hstack((A, B))
# Hz = np.hstack((B.T, A.T))

# # Couple of tests:
# test_stabs_commute(Hx, Hz)
