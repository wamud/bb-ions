''' noisefuncs
These are functions that accept a time and return a value of p.
For example p_idle(T = 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience an error with probability 1e-6
'''

import numpy as np

''' p_idle
If a qubit is idling for time t, this function returns what the probability of an error occurring on it will be. This probability can be fed into other functions to say what the error will actually be. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50 seconds'''
def p_idle(t):
    T2 = 50
    p = 0.5 * (1 - np.exp(-t / T2))
    return p


''' p_merge
Depends on length of time merging the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.
'''
def p_merge(t):
    return t


''' p_shuttle
Depends on length of time being shuttled for. Todo: fill this function. At the moment is just a dummy function returning t.
'''
def p_shuttle(t):
    return t

''' p_split
Depends on length of time splitting the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.
'''
def p_split(t):
    return t
