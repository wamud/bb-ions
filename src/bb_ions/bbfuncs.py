'''A Bicycle Bivariate (BB) code [2308.07915] is a CSS code and thus has separate Hx and Hz parity check matrices. 
These functions are for constructing and verifying the parity check matrices and logical operators of a BB code.'''


import numpy as np



'''make_s
makes a cyclic matrix S. This is an idenity matrix with every 1 cyclically shifted to the right by one'''
def make_s(dim):
  s = np.zeros((dim,dim), dtype = int)

  for i in range(dim):
    s[i % dim, (i + 1) % dim] = 1

  return s

'''make_x
Makes the matrix x, which is x := S_l ⊗ I_m'''
def make_x(l, m):

  s_l = make_s(l)
  ident_m = np.eye(m, dtype = int)

  x = np.kron(s_l, ident_m)
  return x

'''make_y
Makes the matrix y, which is y := I_l ⊗ S_m'''
def make_y(l,m):

  ident_l = np.eye(l, dtype = int)
  s_m = make_s(m)

  y = np.kron(ident_l, s_m)
  return y

'''make_z
Makes the matrix z, which is x * y or, equivalently, S_l ⊗ S_m'''
def make_z(l,m):
  s_l = make_s(l)
  s_m = make_s(m)

  z = np.kron(s_l, s_m)

  return z

'''make_xyz
Uses above functions to make x, y and z'''
def make_xyz(l,m):
  x = make_x(l, m)
  y = make_y(l, m)
  z = make_z(l, m)

  return x, y, z



'''test_stabs_commute
Given X and Z parity check matrices, this function tests that all the stabilisers commute by testing Hx Hz^T = [0] mod 2
(they should have an even number of 1's in corresponding positions so give the zero matrix)'''
def test_stabs_commute(Hx, Hz):
  test = (Hx @ Hz.T) % 2 # should give all-zero matrix zeros
  if np.all(test == 0) != True:
    print("Tragédie, tragédie ! Not all the stabilisers commute")


'''verify_ones
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


''' get_nonzero_indices
Returns a list of the column numbers of the maximum element in each row in a matrix.
When working with matrices that contain all zeros apart from one 1 per row (verified with verify_ones),
this function returns a list which contains the column of that one 1 for each row'''
def get_one_positions(matrix):
  cols = np.argmax(matrix, axis = 1) # axis = 1 means search row by row (as opposed to columns)

  if np.all(cols == 0):
    return None

  return cols


'''make_ones_positions_dict
For each matrix in the matrix_dict, this will make an entry in the dictionary ones_positions_dict.
Each entry contains a list of the column numbers (positions) that contains the 1 for each row.
E.g. if fed A1, then ones_positions_dict['A1'][0] is the position of the 1 in the zeroeth row of A1
'''
def make_ones_positions_dict(A1, A2, A3, B1, B2, B3):

  matrix_dict = {'A1': A1, 'A1T': A1.T, 'A2': A2, 'A2T': A2.T, 'A3': A3, 'A3T': A3.T,
            'B1': B1, 'B1T': B1.T, 'B2': B2, 'B2T': B2.T, 'B3': B3, 'B3T': B3.T }

  ones_positions_dict = {}
  for name, mat in matrix_dict.items():
    ones_positions_dict[name] = get_one_positions(mat)

  return ones_positions_dict


# Using autqec

def find_logical_ops_and_assert_anticommute(n, k, Hx, Hz):

  zeros = np.zeros_like(Hx)
  H_symp = np.array(np.vstack((np.hstack((Hx,zeros)),np.hstack((zeros,Hz)))),dtype=int)

  H_symp_rref, _, transform_rows, transform_cols = rref_mod2(H_symp)
  H_symp_rref = H_symp_rref[~np.all(H_symp_rref == 0, axis=1)]
  H_symp_rref_og_basis = H_symp_rref@inv_mod2(transform_cols)
  assert H_symp_rref_og_basis.shape[0] == n-k
  assert H_symp_rref_og_basis.shape[1] == 2*n

  G, LX_symplectic, LZ_symplectic, D = compute_standard_form(H_symp_rref_og_basis)

  # We have a CSS code so cut off the Z and X parts of the symplectic representations of Lx and Lz:
  Lx = LX_symplectic[:, :n]
  Lz = LZ_symplectic[:, n:]

  # Verify anticommutation of logical operators:
  anticommute_matrix = (Lx @ Lz.T) % 2
  ident_matrix = np.eye(Lx.shape[0])

  # Assert that the anticommute matrix is the identity:
  assert(np.array_equal(anticommute_matrix, ident_matrix)) 
  ## Actually doesn't need to be exactly identity to mean they all anticommute (up to being
  #  multiplied by eachother), a binary matrix of rank 4 is equivalent.
  # However I'd guess there's some benefit to it being the identity.

  return Lx, Lz


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