import unittest
import numpy as np
from bb_ions.bb_code import bb_check_matrices, generate_matrix_polynomial
#from bb_ions.bbparamfuncs import get_code_params

class BBCodeTest(unittest.TestCase):
    """
    Tests for bb code generation
    """
    def test_gen_mat_poly(
        self,
        x=None,
        y=None 
    ):
        np.random.seed(1)

        if x is None:
            x = np.random.randint(1000, size=(10, 10))
        if y is None:
            y = np.random.randint(1000, size=(10, 10))
        
        pows=((0, 2, 3), (2, 1, 4))
        
        p = generate_matrix_polynomial(x, y, pows)
        
        p_targ = y @ y + x @ x @ y + x @ x @ x @ y @ y @ y @ y
        
        assert np.array_equal(p, p_targ)
    
    def test_bb_check_matrixes(
        self,
        l=12,
        m=12,
        left_pow=((3, 0, 0), (0, 2, 7)),
        right_pow=((0, 1, 2), (3, 0, 0))
    ):
        hx, hz = bb_check_matrices(l, m, left_pow, right_pow)

       # aij = tuple(zip(*left_pow))
       # bij = tuple(zip(*right_pow))

       # code = get_code_params(l, m, aij, bij)

       # assert np.array_equal(code.Hx, hx)
       # assert np.array_equal(code.Hz, hz)

if __name__ == "__main__":
    unittest.main()
