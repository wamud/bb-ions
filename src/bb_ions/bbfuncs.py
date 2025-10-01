''' bbfuncs
A Bicycle Bivariate (BB) code [2308.07915] is a CSS code and thus has separate Hx and Hz parity check matrices. 
These functions are for constructing the parity check matrices of a BB code.'''

import numpy as np

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



''' test_stabs_commute
Given X and Z parity check matrices, this function tests that all the stabilisers commute by testing Hx Hz^T = [0] mod 2
(they should have an even number of 1's in corresponding positions so give the zero matrix)'''
def test_stabs_commute(Hx, Hz):
  test = (Hx @ Hz.T) % 2 # should give all-zero matrix zeros
  if np.all(test == 0) != True:
    print("Tragédie, tragédie ! Not all the stabilisers commute")


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
# verify_ones(A1, A3, B1, B2, B3)