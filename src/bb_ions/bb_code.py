"""
This module defines the BBCode class which is used to generate check matrices 
and logical operators for a bivariate-bycycle code.

Ref: 
"""

from typing import Iterable, Tuple
import numpy as np
import json
from numpy.linalg import matrix_power
from bposd.css import css_code
from bposd.css_decode_sim import css_decode_sim

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
        "error_rate": 0.1
    }

    def __init__(
        self,
        l: int,
        m: int,
        left_pow: Iterable[Iterable[int]],
        right_pow: Iterable[Iterable[int]],
        estimate_distance = False
    ):
        """
        args
        ----
        l: 
        m: 
        left_pow:
        right_pow:
        """
        self.l = l
        self.m = m
        self.left_pow = left_pow
        self.right_pow = right_pow
        self.Hx, self.Hz = self.generate_check_mat()
        self.N, self.K = self.generate_nk()
        self.distance = None
        if estimate_distance:
            self.distance = self.estimate_distance()

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
        s_l = np.roll(np.eye(self.l, dtype=int), 1,  axis=1)
        
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
        code = css_code(hx=self.Hx, hz=self.Hz)
        return code.N, code.K

    def estimate_distance(self):
        """
        Estimate distance of the code
        """
        
        if self.distance is not None:
            return self.distance

        bb5 = css_decode_sim(
            hx=self.Hx, hz=self.Hz, **self.OSD_OPTIONS
        )

        results = json.loads(bb5.output_dict())

        dmax = results["min_logical_weight"]
        
        self.distance = dmax
        
        return dmax

    def generate_logical_ops(self):
        pass

    def __str__(self):
        A = ' + '.join([f"x^{i} * y^{j}" for i, j in zip(*self.left_pow)])
        B = ' + '.join([f"x^{i} * y^{j}" for i, j in zip(*self.right_pow)])
        n = self.N
        k = self.K

        if self.distance is None:
            d = '?'
        else:
            d = self.distance
        
        return '\n'.join(
            [f"[[{n}, {k}, <={d}]]",
             f"l = {self.l}",
             f"m = {self.m}",
             f"A = {A}",
             f"B = {B}"]
        )
    
