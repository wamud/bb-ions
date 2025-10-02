''' noisefuncs
Hardware-inspired noise and error functions relevant to our architecture
'''

import numpy as np

''' p_idle
If a qubit is idling for time t, this function returns what the probability of an error occurring on it will be. This probability can be fed into other functions to say what the error will actually be. For example p_idle(T = 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience an error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50 seconds'''
def p_idle(t):
    T2 = 50
    p = 0.5 * (1 - np.exp(-t / T2))
    return p


''' p_merge
Depends on length of time merging the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.
'''
def p_merge(t):
    return t



''' p_shift
Depends on length of time being shuttled for. Todo: fill this function properly. At the moment will just say it's proportional to the difference between the index of the previous leg the modules were in and the next.
'''
def p_shift(t_shift):
    
    p = 1e-3 * t_shift
    
    return p

''' apply_shift_error
Accepts the stim circuit to have operations appended to, the register of qubits that require a cyclic shift error and the probability of a cyclic shift error.
Appends depolarising noise onto the qubits being cyclically shifted with strength p'''
def apply_shift_error(circuit, register, p):
    if p != 0:
        circuit.append("DEPOLARIZE1", register, p)


''' p_split
Depends on length of time splitting the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.
'''
def p_split(t):
    return t

