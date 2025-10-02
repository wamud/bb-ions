''' modulefuncs
New functions since making modular according to Tham, Ye, Khait, Gamble, Delfosse [2508.01879]: Quantum memories over a 2 × L array of qubit modules'''

from .kfuncs import *

''' findIJ
Given a bivariate polynomial, for example 
P = x^0 + y + x^2*y^2  
And writing tuples (i, j) for each term x^i⋅y^j, for example
P = [(0, 0), (0, 1), (2, 2)]
This function returns
I(P) = list of powers of x in P
J(P) = list of powers of y in P
'''
def findIJ(P):
    IP = [P[k][0] for k in range(len(P))]
    JP = [P[k][1] for k in range(len(P))]
    return IP, JP


''' new_make_registers
Makes lists of qubit indices dividing n = 2lm qubits evenly into qA, qB, qC, qD.
  qA: qubits 0 to n/2 - 1
  qB: qubits n/2 to n - 1
  qC: qubits n to 3n/2 - 1 
  qD: qubits 3n/2 to 2n - 1
This is the same as the make_registers function, it's just that it uses the convtok within it'''
def new_make_registers(l, m):
  qA = [convtok(l, m, 0, v, w) for v in range(l) for w in range(m)]
  qB = [convtok(l, m, 1, v, w) for v in range(l) for w in range(m)]
  qC = [convtok(l, m, 2, v, w) for v in range(l) for w in range(m)]
  qD = [convtok(l, m, 3, v, w) for v in range(l) for w in range(m)]
  return qA, qB, qC, qD
  