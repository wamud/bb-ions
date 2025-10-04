''' noisefuncs
Hardware-inspired noise and error functions relevant to our architecture
'''

import numpy as np

''' p_init
The probability of initalising to an orthogonal eigenstate. Depends on time taken to initalise (t_init). todo: fill with more than dummy function.'''
def p_init(t_init):
    p = t_init
    return p

'''p_had
The probability of an error occuring on the hadamard gate. For what error is actually applied see 'hadamard' in circfuncs.'''
def p_had(t_had):
    p = t_had
    return p


''' p_shuttle 
If shuttling a qubit it will experience an error associated with acceleration, time of shuttle, number of T junctions it passes through and deceleration. 
todo: make this function more than a dummy'''
def p_shuttle(t_shuttle):
    p = t_shuttle
    return p

''' apply_shuttle_error
Applies a depolarising noise channel of strength p to qubits in 'register' and in stim circuit 'circuit'. todo: make more accurate noise model once we have the info'''
def apply_shuttle_error(circuit, register, t_shuttle):
    p = p_shuttle(t_shuttle)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)

''' p_idle
If a qubit is idling for time t, this function returns what the probability of an error occurring on it will be. This probability can be fed into other functions to say what the error will actually be. For example p_idle(T = 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience an error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50 seconds'''
def p_idle(t):
    T2 = 50
    p = 0.5 * (1 - np.exp(-t / T2))
    return p


''' idle
Adds dephasing (Z) noise of strength p (a Z operation is applied with probability p) to qubits in 'register.'''
def idle(circuit, register, t_idle = 0):
  p = p_idle(t_idle)
  if p > 0:
    circuit.append("Z_ERROR", register, p)


''' p_shift
todo: make more than dummy function'''
def p_shift(t_shift):
    p = 1e-3 * t_shift
    return p


''' apply_shift_error
todo: make more than dummy function'''
def apply_shift_error(circuit, register, t_shift):
    p = p_shift(t_shift)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p) 

''' p_merge
Depends on length of time merging the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.'''
def p_merge(t):
    p = t
    return p

''' apply_merge_error
Once check modules have been cyclically shifted to the data qubit module they need to interact with, we simulate merging their coulomb potentials'''
def apply_merge_error(circuit, register, t_merge):
    p = p_merge(t_merge)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)


'''p_cnot
todo: fill '''
def p_cnot(t_cnot):
    p = t_cnot
    return p


'''p_cz
todo: fill '''
def p_cz(t_cz):
    p = t_cz
    return p



''' p_split
todo: make more than dummy function'''
def p_split(t):
    p = t
    return p

''' apply_split_error
Simulating splitting the coulomb potentials of check qubit and data qubit modules by appling an error probability to them. At the moment is just depolarising noise of strength p = t_split'''
def apply_split_error(circuit, register, t_split):
    p = p_split(t_split)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)

''' p_meas
Probability of a measurement error given the time of measurement was t_meas. At the moment just returns p equal to t_meas'''
def p_meas(t_meas):
        p = t_meas
        return p
