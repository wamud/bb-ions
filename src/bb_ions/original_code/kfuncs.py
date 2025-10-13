"""
Here we have functions to convert between indexing 2lm qubits using 2lm
integers k = {0, 1, 2, ... , 2lm - 1} and using a row / column notation to
align with the qubit indexing in [2508.01879] which enables its Algorithm 2 for
the syndrome extraction circuit.

We split the qubits into two groups, u = 0 and u = 1. Each group has 2lm
qubits. We split the 2lm qubits into rows of size m. Having m columns per row
implies we have l rows. We can now refer to a qubit by group / row / column (u,
v, w).

We convert from k = {0, 1, 2, ... , 2lm - 1} to row / column (v / w) notation
by defining (v, w) = (⌊k/m⌋ , k mod m). These are standard conversion
definitions to arrange into groups of m -- note the value of v only goes up by
1 every m elements and w cycles through m elements, hence a row and column.

To convert from (u, v, w) to the original qubit index is k = ulm + vm + w.
"""

import math


def convtorowcol(m, k):
    """ 
    Converts an integer k to row / column notation where each row has m
    elements. This is done using k = (v, w) := (⌊k/m⌋ , k mod m).
    """

    v = math.floor(k / m)
    w = k % m
    return v, w



def convtouvw(l, m, k):
    """ 
    If k indexes qubits such that k ∈ {0, 1, ...} and we want to arrange k into
    groups of lm qubits, where each group is an array with l rows and m
    columns, then this function returns (u, v, w) where 
    - u indicates which group / array of lm qubits it's in (0th, 1st, 2nd etc.
      -- can choose 0th & 1st to be data qubits)
    - v is the row within the array 
    - w is the column within the array
    """

    u = math.floor(k / (l * m))
    this_arrays_k = k - u * l * m
    v, w = convtorowcol(m, this_arrays_k)
    return u, v, w



def conv_vw_to_k(v, w, m):
    """ 
    This finds the index (in reading order) of an element at row / column (v,
    w) in an array where each row has m elements
    """

    assert w < m, "w must be less than m -- there are only m columns per row"
    k = v * m + w
    return k



def convtok(l, m, u, v, w):
    """ 
    Converts from the tuple (u, v, w), where u indicates which block of lm
    qubits a qubit is in and v and w are its row and column within that block,
    to an index k, which simply goes from 0 up to number of qubits in reading
    order
    """

    assert w < m, "w must be less than m -- there are only m columns per row"
    assert v < l, "v must be less than l -- there are only l rows per group"
    k = u * l * m + v * m + w
    return k
