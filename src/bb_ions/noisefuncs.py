''' noisefuncs
Hardware-inspired noise and error functions relevant to our architecture'''

import numpy as np


# class Error:
#     def __init__(self, p, type):
#         self.p = p
#         self.type = type

# class Errors:
#     def __init__(self, p, RZ, RX):
#         self.RZ = Error(p, 'X_ERROR')
#         self.RX = Error(p, 'Z_ERROR')

class Errors:
    def __init__(self, p, 

                p_RZ = None, RZ_error = 'X_ERROR',
                p_RX = None, RX_error = 'Z_ERROR',
                p_H = None, H_error = 'DEPOLARIZE1',
                p_MZ = None,
                p_MX = None,

                p_CNOT = None, CNOT_error = 'DEPOLARIZE2',
                p_CZ = None, CZ_error = 'DEPOLARIZE2',

                p_shift_const = None, shift_error = 'DEPOLARIZE1',
                p_shuttle = 0, shuttle_error = 'DEPOLARIZE1',
                p_merge = 0, merge_error = 'DEPOLARIZE1',
                p_split = 0, split_error = 'DEPOLARIZE1'

                ):
        
        # Single-qubit operations:
        
        self.p_RZ = p_RZ or p
        self.RZ_error = RZ_error

        self.p_RX = p_RX or p
        self.RX_error = RX_error 

        self.p_H = p_H or p
        self.H_error = H_error 

        self.p_MZ = p_MZ or p

        self.p_MX = p_MX or p

        # Two-qubit operations:
        
        self.p_CNOT = p_CNOT or p
        self.CNOT_error = CNOT_error

        self.p_CZ = p_CZ or p
        self.CZ_error = CZ_error

        # Additional operations:

        self.p_shift_const = p_shift_const or p # setting p_shift_const to p, though note within apply_shift_error this is multiplied by the length of the cyclic shift hence it's 'p_shift_const'
        self.shift_error = shift_error

        # Default values of additional operations is zero:
        self.p_shuttle = p_shuttle
        self.shuttle_error = shuttle_error

        self.p_merge = p_merge 
        self.merge_error = merge_error
        
        self.p_split = p_split 
        self.split_error = split_error




class Idlings:
    def __init__(self, p, 
                RZ = None, 
                RX = None, 
                H = None, 
                MZ = None, 
                MX = None,
                CNOT = None,
                CZ = None, 
                ):
        
        # Idling during single-qubit operations:
        
        self.RZ = [RZ[0] or p, RZ[1] or 'DEPOLARIZE1']
        self.RZ = [ RZ[0] or p, RZ[1] or 'DEPOLARIZE1']
        self.RX = [ RX[0] or p, RX[1] or 'DEPOLARIZE1']
        self.H = [ H[0] or p, H[1] or 'DEPOLARIZE1']
        self.MZ = [ MZ[0] or p, MZ[1] or 'DEPOLARIZE1']
        self.MX = [ MX[0] or p, MX[1] or 'DEPOLARIZE1']

        # Idling during two-qubit operations:
        self.CNOT = [ CNOT[0] or p, CNOT[1] or 'DEPOLARIZE1']
        self.CZ = [ CZ[0] or p, CZ[1] or 'DEPOLARIZE1']




class NoiseTimes:
    def __init__(self, t_init, t_had, t_merge, t_split, t_cnot, t_cz, t_shuttle, t_shift_const, t_meas, t_idle, t_idle_meas):
        self.t_init = t_init
        self.t_had = t_had
        self.t_merge = t_merge
        self.t_split = t_split
        self.t_cnot = t_cnot
        self.t_cz = t_cz
        self.t_shuttle = t_shuttle
        self.t_shift_const = t_shift_const
        self.t_meas = t_meas
        self.t_idle = t_idle
        self.t_idle_meas = t_idle_meas




''' make_uniform_noisetimes
Create an object noisetimes of class NoiseTimes to store the lengths of time each operation takes so they can be used when applying noise. This sets them all uniformly to the same input to this function t. The object noisetimes can be modified afterwards to set specific times'''
def make_uniform_noisetimes(t):
    t_init = t
    t_had = t
    t_merge = t
    t_split = t
    t_cnot = t
    t_cz = t
    t_shuttle = t
    t_shift_const = t
    t_meas = t
    t_idle = t
    t_idle_meas = t

    noisetimes = NoiseTimes(t_init, t_had, t_merge, t_split, t_cnot, t_cz, t_shuttle, t_shift_const, t_meas, t_idle, t_idle_meas)
    
    return noisetimes


''' make_longchain_noisetimes
Create an object noisetimes of class NoiseTimes to store the lengths of time each operation takes so they can be used when applying noise. Sets according to Ye Delfosse longchain: 2503.2207'''
def make_longchain_noisetimes(t):
    
    t_init = t / 10
    t_had = t / 10
    t_cnot = t
    t_cz = t
    t_meas = t / 10
    
    t_idle_meas = 30 * t / 100
    t_idle = t / 100

    # Our operations:
    t_shuttle = t / 10
    t_shift_const = t / 10
    t_merge = t / 10
    t_split = t / 10

    noisetimes = NoiseTimes(t_init, t_had, t_merge, t_split, t_cnot, t_cz, t_shuttle, t_shift_const, t_meas, t_idle, t_idle_meas)
    
    return noisetimes


''' make_uniform_agnostic_noisetimes
Sets noise of all operations to the input t except removes the noise introduced by our hardware proposal, namely merge, split, shuttle, shift.'''
def make_uniform_agnostic_noisetimes(t):
    t_init = t
    t_had = t
    t_cnot = t
    t_cz = t
    t_meas = t
    t_idle = t
    t_idle_meas = t

    # Our hardware proposal:
    t_merge = 0
    t_split = 0
    t_shuttle = 0
    t_shift_const = 0

    noisetimes = NoiseTimes(t_init, t_had, t_merge, t_split, t_cnot, t_cz, t_shuttle, t_shift_const, t_meas, t_idle, t_idle_meas)
    
    return noisetimes


''' make_longchain_agnostic_noisetimes
Sets noise of all operations as per Ye Delfosse longchain paper 2503.2207. Also removes noise introduced by our hardware proposal, namely merge, split, shuttle, shift are set to zero.'''
def make_longchain_agnostic_noisetimes(t):
    t_init = t / 10
    t_had = t / 10
    t_cnot = t
    t_cz = t
    t_meas = t / 10
    
    t_idle_meas = 30 * t / 100
    t_idle = t / 100

    # Our hardware proposal:
    t_merge = 0
    t_split = 0
    t_shuttle = 0
    t_shift_const = 0

    noisetimes = NoiseTimes(t_init, t_had, t_merge, t_split, t_cnot, t_cz, t_shuttle, t_shift_const, t_meas, t_idle, t_idle_meas)
    
    return noisetimes



''' p_init
The probability of initalising to an orthogonal eigenstate. Depends on time taken to initalise (t_init). todo: fill with more than dummy function.'''
def p_init(t):
    p = t
    return p

'''p_had
The probability of an error occuring on the hadamard gate. For what error is actually applied see 'hadamard' in circfuncs.'''
def p_had(t):
    p = t
    return p

''' p_shuttle 
If shuttling a qubit it will experience an error associated with acceleration, time of shuttle, number of T junctions it passes h e e e  through and deceleration. 
todo: make this function more than a dummy'''
def p_shuttle(t_shuttle):
    p = t_shuttle
    return p

''' p_shift
todo: make more than dummy function'''
def p_shift(t_shift):
    p = t_shift
    return p

''' p_merge
Depends on length of time merging the coulomb potential of two ion modules takes. Todo: fill this function. At the moment is just a dummy function returning t.'''
def p_merge(t):
    p = t
    return p

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

''' p_meas
Probability of a measurement error given the time of measurement was t_meas. At the moment just returns p equal to t_meas'''
def p_meas(t_meas):
        p = t_meas
        return p

''' p_idle'''
def p_idle(t_idle):
    p = t_idle
    return p 

# ''' p_idle
# If a qubit is idling for time t, this function returns what the probability of an error occurring on it will be. This probability can be fed into other functions to say what the error will actually be. For example p_idle(T = 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience an error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50 seconds'''
# def p_idle(t):
#     T2 = 50
#     p = 0.5 * (1 - np.exp(-t / T2))
#     return p



''' apply_shuttle_error
Applies a depolarising noise channel of strength p to qubits in 'register' and in stim circuit 'circuit'. todo: make more accurate noise model once we have the info'''
def apply_shuttle_error(circuit, register, t_shuttle):
    p = p_shuttle(t_shuttle)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)


# ''' idle
# Adds dephasing (Z) noise of strength p (a Z operation is applied with probability p) to qubits in 'register.'''
# def idle(circuit, register, t_idle = 0):
#   p = p_idle(t_idle)
#   if p > 0:
#     circuit.append("Z_ERROR", register, p)


''' idle
Adds depolarising noise of strength p_idle(t_idle) to qubits in register'''
def idle(circuit, register, t_idle = 0):
  p = p_idle(t_idle)
  if p > 0:
    circuit.append("DEPOLARIZE1", register, p)

''' apply_shift_error
todo: make more than dummy function'''
def apply_shift_error(circuit, register, t_shift):
    p = p_shift(t_shift)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p) 


''' apply_merge_error
Once check modules have been cyclically shifted to the data qubit module they need to interact with, we simulate merging their coulomb potentials'''
def apply_merge_error(circuit, register, t_merge):
    p = p_merge(t_merge)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)

''' apply_split_error
Simulating splitting the coulomb potentials of check qubit and data qubit modules by appling an error probability to them. At the moment is just depolarising noise of strength p = t_split'''
def apply_split_error(circuit, register, t_split):
    p = p_split(t_split)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)
