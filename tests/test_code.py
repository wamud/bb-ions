import unittest
import numpy as np
from bb_ions.code import generate_matrix_polynomial, BBCode
# from bb_ions.bbparamfuncs import get_code_params


class BBCodeTest(unittest.TestCase):
    """
    Tests for bb code generation
    """

    def test_gen_mat_poly(self, x=None, y=None):
        np.random.seed(1)

        if x is None:
            x = np.random.randint(1000, size=(10, 10))
        if y is None:
            y = np.random.randint(1000, size=(10, 10))

        pows = ((0, 2, 3), (2, 1, 4))

        p = generate_matrix_polynomial(x, y, pows)

        p_targ = y @ y + x @ x @ y + x @ x @ x @ y @ y @ y @ y

        assert np.array_equal(p, p_targ)

    def test_generate_check_mat(
        self,
        l=12,
        m=12,
        left_pow=((3, 0, 0), (0, 2, 7)),
        right_pow=((0, 1, 2), (3, 0, 0)),
    ):
        bb = BBCode(l, m, left_pow, right_pow)
        hx, hz = bb.generate_check_mat()

        # aij = tuple(zip(*left_pow))
        # bij = tuple(zip(*right_pow))

        # code = get_code_params(l, m, aij, bij)

        # assert np.array_equal(code.Hx, hx)
        # assert np.array_equal(code.Hz, hz)

    def test_generate_nk(
        self,
        l=12,
        m=12,
        left_pow=((3, 0, 0), (0, 2, 7)),
        right_pow=((0, 1, 2), (3, 0, 0)),
        targ_n=288,
        targ_k=12,
    ):
        bb = BBCode(l, m, left_pow, right_pow)

        n, k = bb.generate_nk()

        assert n == targ_n
        assert k == targ_k

    def test_estimate_distance(
        self,
        l=12,
        m=12,
        left_pow=((3, 0, 0), (0, 2, 7)),
        right_pow=((0, 1, 2), (3, 0, 0)),
    ):
        bb = BBCode(l, m, left_pow, right_pow)
        print(bb.estimate_distance())


if __name__ == "__main__":
    unittest.main()
