'''Here we have functions to convert between indexing 2lm qubits
using 2lm integers k = {0, 1, 2, ... , 2lm - 1} and using a 
row / column notation, organising the qubits into two groups 
of l rows and m columns such that we can refer to qubit k by 
group / row / column (u, v, w).

This can be done by defining (v, w) = (⌊k/m⌋ , k mod m).
(These are standard conversion definitions to arrange into groups of m --
note the value of v only goes up by 1 every m elements and w cycles through m
elements, hence a row and column)

To convert back is k = ulm + vm + w

We are doing this to align with the qubit indexing in [2508.01879] and enable its
Algorithm 2 for the syndrome extraction circuits of a BB code arranged
into modules, where each module has qubits all with the same column number w.
'''

import math

'''convtorowcol
Converts an integer k to row / column notation where each row has m elements.
This is done using k = (v, w) := (⌊k/m⌋ , k mod m).'''
def convtorowcol(m, k):
    v = math.floor(k/m)
    w = k % m
    return v, w


'''convtok
This finds the index (in reading order) of an element at row / column (v, w) in an array
where each row has m elements'''
def convtok(v, w, m):
    k = v * m + w
    return k

'''convtouvw
If k indexes qubits such that k ∈ {0, 1, ...} and we want to arrange k into groups
of lm qubits, where each group is an array with l rows and m columns, then this function returns (u, v, w) where 
    u indicates which group / array of lm qubits its in
    v is the row within the array
    w is the column within the array'''
def convtouvw(l, m, k):
    u = math.floor(k / (l * m) )
    this_arrays_k = k - u * l * m 
    v, w = convtorowcol(m, this_arrays_k)
    return u, v, w

'''conv_uvw_to_k
Converts from the tuple (u, v, w), where u indicates which block of lm qubits a qubit is in and v and w are its row and column within that block, to an index k, which simply goes from 0 up to number of qubits in reading order'''
def conv_uvw_to_k(l, m, u, v, w):
    k = u * l * m   +   v * m   +   w
    return k

