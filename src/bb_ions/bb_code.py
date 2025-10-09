"""
This module defines the BBCode class which is used to generate check matrices
and logical operators for a bivariate-bycycle code.

Ref:
"""

from typing import Iterable, Tuple
import numpy as np
import json
from functools import cached_property
from numpy.linalg import matrix_power
from bposd.css import css_code
from bposd.css_decode_sim import css_decode_sim
from autqec.utils.qec import rref_mod2, inv_mod2, compute_standard_form


def generate_matrix_polynomial(
    x: np.ndarray, y: np.ndarray, poly_pow: Iterable[Iterable[int]]
) -> np.ndarray:
    """
    Generate a matrix polynomial with variables x and y.

    Assumes x and y are square matrices
    """
    return sum(matrix_power(x, i) @ matrix_power(y, j) for i, j in zip(*poly_pow))


class CodeQubits:
    def __init__(self, n, k):
        self.physical = n
        self.logical = k


class CodeOperators:
    def __init__(self, x, z):
        self.X = x
        self.Z = z


class BBCode:
    """
    Defines a bivariate-bycycle code
    """

    OSD_OPTIONS = {
        "target_runs": 2000,
        "xyz_error_bias": [1, 1, 1],
        "bp_method": "minimum_sum",
        "ms_scaling_factor": 0.05,
        "osd_method": "osd_cs",
        "osd_order": 4,
        "channel_update": None,
        "seed": 42,
        "max_iter": 9,
        "error_bar_precision_cutoff": 1e-6,
        "tqdm_disable": True,
        "error_rate": 0.1,
    }

    def __init__(
        self,
        l: int,
        m: int,
        left_pow: Iterable[Iterable[int]],
        right_pow: Iterable[Iterable[int]],
        estimate_distance=False,
    ):
        """
        args
        ----
        l:
        m:
        left_pow: Iterable of x and y powers for left matrix polynomial.
                  x^ay^b + x^cy^d + x^ey^f is represented by
                  `left_pow` = ((a, c, e), (b, d, f))
        right_pow: Iterable of x and y powers for right matrix polynomial
                  x^ay^b + x^cy^d + x^ey^f is represented by
                  `right_pow` = ((a, c, e), (b, d, f))
        """
        self.l = l
        self.m = m
        self.left_pow = tuple(left_pow)
        self.right_pow = tuple(right_pow)

        self.distance = None
        if estimate_distance:
            self.distance = self.estimate_distance()

    @cached_property
    def qubits(self):
        n, k = self.generate_nk()
        return CodeQubits(n, k)

    @cached_property
    def check_operators(self):
        Hx, Hz = self.generate_check_mat()
        return CodeOperators(Hx, Hz)

    @cached_property
    def logical_operators(self):
        Lx, Lz = self.generate_logical_ops()

        return CodeOperators(Lx, Lz)

    def generate_check_mat(self):
        """
        Generates Hx and Hz for the BB code

        args

        Returns
        -------
        (Hx, Hz) the X and Z parity check matrices
        ----
        """
        # Identity matrix cyclically shifted by 1 column
        s_l = np.roll(np.eye(self.l, dtype=int), 1, axis=1)

        # generate variables
        x = np.kron(s_l, np.eye(self.m, dtype=int))
        y = np.kron(np.eye(self.m, dtype=int), s_l)

        # generate polynomial
        A = generate_matrix_polynomial(x, y, self.left_pow)
        B = generate_matrix_polynomial(x, y, self.right_pow)

        Hx = np.hstack((A, B))
        Hz = np.hstack((B.T, A.T))

        # check that all stabilizer checks commute
        assert not np.any((Hx @ Hz.T) % 2)

        return Hx, Hz

    def generate_nk(self):
        code = css_code(hx=self.check_operators.X, hz=self.check_operators.Z)
        return code.N, code.K

    def estimate_distance(self):
        """
        Estimate distance of the code
        """

        if self.distance is not None:
            return self.distance

        bb5 = css_decode_sim(
            hx=self.check_operators.X, hz=self.check_operators.Z, **self.OSD_OPTIONS
        )

        results = json.loads(bb5.output_dict())

        dmax = results["min_logical_weight"]

        self.distance = dmax

        return dmax

    def generate_logical_ops(self):
        """
        Generate the logical operators of the code.
        """

        # Convert to symplectic form
        zeros = np.zeros_like(self.check_operators.X)
        H_symp = np.array(
            np.vstack(
                (
                    np.hstack((self.check_operators.X, zeros)),
                    np.hstack((zeros, self.check_operators.Z)),
                )
            ),
            dtype=int,
        )

        # Row reduce and find logical operators
        H_symp_rref, _, transform_rows, transform_cols = rref_mod2(H_symp)
        H_symp_rref = H_symp_rref[~np.all(H_symp_rref == 0, axis=1)]
        H_symp_rref_og_basis = H_symp_rref @ inv_mod2(transform_cols)
        assert (
            H_symp_rref_og_basis.shape[0] == self.qubits.physical - self.qubits.logical
        )
        assert H_symp_rref_og_basis.shape[1] == 2 * self.qubits.physical

        G, LX_symplectic, LZ_symplectic, D = compute_standard_form(H_symp_rref_og_basis)

        # We have a CSS code so cut off the empty Z and X parts of the symplectic
        # representations of Lx and Lz:
        Lx = LX_symplectic[:, : self.qubits.physical]
        Lz = LZ_symplectic[:, self.qubits.physical :]

        # Verify anticommutation of logical operators and that they're canonical:
        assert np.array_equal((Lx @ Lz.T) % 2, np.eye(Lx.shape[0]))

        return Lx, Lz

    def __str__(self):
        A = " + ".join([f"x^{i} * y^{j}" for i, j in zip(*self.left_pow)])
        B = " + ".join([f"x^{i} * y^{j}" for i, j in zip(*self.right_pow)])
        n = self.qubits.physical
        k = self.qubits.logical

        if self.distance is None:
            d = " "
        else:
            d = "<=" + str(self.distance)

        return "\n".join(
            [
                f"[[{n}, {k}, {d}]]",
                f"l = {self.l}",
                f"m = {self.m}",
                f"A = {A}",
                f"B = {B}",
            ]
        )
