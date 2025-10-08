"""
New functions since making modular according to Tham, Ye, Khait, Gamble,
Delfosse [2508.01879]: Quantum memories over a 2 × L array of qubit modules
"""

from .kfuncs import *


def findIJ(P):
    """ 
    Given a bivariate polynomial, for example P = x^0 + y + x^2*y^2  And
    writing tuples (i, j) for each term x^i⋅y^j, for example P = [(0, 0), (0,
    1), (2, 2)] This function returns I(P) = list of powers of x in P J(P) =
    list of powers of y in P 
    """

    IP = [P[k][0] for k in range(len(P))]
    JP = [P[k][1] for k in range(len(P))]
    return IP, JP
