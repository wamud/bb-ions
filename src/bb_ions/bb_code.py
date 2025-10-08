"""
This module generates the parity check matrices for a bivariate-bycycle code

Ref: 
"""

from typing import Iterable, Tuple
import numpy as np
from numpy.linalg import matrix_power

def generate_matrix_polynomial(
        x: np.ndarray,
        y: np.ndarray,
        poly_pow: Iterable[Iterable[int]]
) -> np.ndarray:
    """
    Generate a matrix polynomial with variables x and y.

    Assumes x and y are square matrices
    """
    return sum(
        matrix_power(x, i) @ matrix_power(y, j)
        for i, j in zip(*poly_pow)
    )

def bb_check_matrices(
        l: int,
        m: int,
        left_powers: Iterable[Iterable[int]],
        right_powers: Iterable[Iterable[int]]) -> Tuple[np.ndarray]:
    """
    Generates Hx and Hz for the BB code

    args

    Returns
    -------
    (Hx, Hz) the X and Z parity check matrices 
    ----
    """
    # Identity matrix cyclically shifted by 1 column
    s_l = np.roll(np.eye(l, dtype=int), 1,  axis=1)
    
    # generate variables
    x = np.kron(s_l, np.eye(m, dtype=int))
    y = np.kron(np.eye(m, dtype=int), s_l)
    
    # generate polynomial
    A = generate_matrix_polynomial(x, y, left_powers) 
    B = generate_matrix_polynomial(x, y, right_powers) 
    
    Hx = np.hstack((A, B))
    Hz = np.hstack((B.T, A.T))
    
    # check that all stabilizer checks commute
    assert not np.any((Hx @ Hz.T) % 2)
    
    return Hx, Hz 
