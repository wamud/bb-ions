''' modulefuncs
We will split 2lm data qubits and 2lm check qubits into m 'modules', each of size 2l. This will be done by first organsing 2lm qubits into two groups (left and right data qubits), each with l rows of size m and referring to each qubit by its row and column, (v, w). This was done using kfuncs. We now group qubits with the same column index w into the same module.'''


# Actually not sure I need these functions....


# from kfuncs import *

# ''' make_data_modules
# In a BB code, the first lm data qubits are "left" data qubits (acted on by A in the left side of Hx) and the second lm data qubits are the "right" data qubits (acted on by B in the right side of Hx).
# This uses kfuncs to go through each qubit and group them according to their (u, v, w) values into modules

# We will say from qubits 0 to 2lm - 1 are data qubits. From 0 to lm - 1 are left data qubits, and from lm to 2lm - 1 are right data qubits.
# '''
# def make_data_modules(l, m):
#     for k in range(2 * l * m):

