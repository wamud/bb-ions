''' noisefuncs'''

import numpy as np


class Error:
    def __init__(self, operation, p):
        self.op = operation # the error operation
        self.p = p # its probability

def default_errors(p = 0.001):
    """
    Defines default gates and error rates.
    """
    errors = {

        "RZ" : Error("DEPOLARIZE1", p / 10), # note is depolarize (as per longcahin paper) as opposed to X_ERROR

        "RX" : Error("DEPOLARIZE1", p / 10), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", p / 10),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p / 10),
        "MX" : Error("Z_ERROR", p / 10),

        "shuttle" : Error("DEPOLARIZE1", p / 10),
        "merge" : Error("DEPOLARIZE1", p / 10),
        "split" : Error("DEPOLARIZE1", p / 10),
        "shift" : Error("DEPOLARIZE1", p / 10),

        "shift_constant" : p / 10, 
        
    }
    
    return errors

def zero_errors():
    """
    Defines zero error rates.
    """
    errors = {

        "RZ" : Error("DEPOLARIZE1", 0), # note is depolarize (as per longcahin paper) as opposed to X_ERROR

        "RX" : Error("DEPOLARIZE1", 0), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", 0),
        "CNOT" : Error("DEPOLARIZE2", 0),
        "CZ" : Error("DEPOLARIZE2", 0),
        "MZ" : Error("X_ERROR", 0),
        "MX" : Error("Z_ERROR", 0),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_constant" : 0, 
        
    }
    
    return errors



def zero_idle_errors(p = 0.001):
    """
    Defines default idling gates and error rates, i.e. the operation and probability of that operation applied to qubits that are idling in a timestep that other qubits are experiencing the key operation.
    """
    idle_during = {
        "RZ" : Error("DEPOLARIZE1", 0),
        "RX" : Error("DEPOLARIZE1", 0), 
        "H" : Error("DEPOLARIZE1", 0),
        "CNOT" : Error("DEPOLARIZE1", 0),
        "CZ" : Error("DEPOLARIZE1", 0),
        "MZ" : Error("DEPOLARIZE1", 0),
        "MX" : Error("DEPOLARIZE1", 0),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_constant" : 0,
    }

    return idle_during

def default_idle_errors(p = 0.001):
    """
    Defines default idling gates and error rates, i.e. the operation and probability of that operation applied to qubits that are idling in a timestep that other qubits are experiencing the key operation.
    """
    idle_during = {
        "RZ" : Error("DEPOLARIZE1", p / 100),
        "RX" : Error("DEPOLARIZE1", p / 100), 
        "H" : Error("DEPOLARIZE1", p / 100),
        "CNOT" : Error("DEPOLARIZE1", p / 100),
        "CZ" : Error("DEPOLARIZE1", p / 100),
        "MZ" : Error("DEPOLARIZE1", 30 * p / 100),
        "MX" : Error("DEPOLARIZE1", 30 * p / 100),

        "shuttle" : Error("DEPOLARIZE1", p / 100),
        "merge" : Error("DEPOLARIZE1", p / 100),
        "split" : Error("DEPOLARIZE1", p / 100),
        "shift" : Error("DEPOLARIZE1", p / 100),

        "shift_constant" : p / 100,
    }

    return idle_during



# ''' p_idle
# If a qubit is idling for time t, this function returns what the probability of an error occurring on it will be. This probability can be fed into other functions to say what the error will actually be. For example p_idle(T = 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience an error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50 seconds'''
# def p_idle(t):
#     T2 = 50
#     p = 0.5 * (1 - np.exp(-t / T2))
#     return p



''' idle
Adds an idling error to qubits in register. The idling error is of class Error.
E.g. idle(circuit, [0, 1], idle_during['MZ']) '''
def idle(circuit, register, error: Error):

  p = error.p
  
  if p > 0:
    circuit.append(error.op, register, p)


''' apply_shuttle_error
Applies a depolarising noise channel of strength p to qubits in 'register' and in stim circuit 'circuit'. todo: make more accurate noise model once we have the info'''
def apply_shuttle_error(circuit, register, errors: dict):

    p = errors['shuttle'].p

    if p > 0:

        circuit.append(errors['shuttle'].op, register, p)



''' apply_shift_error
Applies an error to qubits to simulate them undergoing the cyclic shift required to align check and data qubit modules. Contained in errors is the shift constant. The actual value of p can be fed in to represent longer or shorter cyclic shifts'''
def apply_shift_error(circuit, register, errors):
    
    p = errors['shift'].p
    if p > 0:
        circuit.append(errors['shift'].op, register, p) 


''' apply_merge_error
Once check modules have been cyclically shifted to the data qubit module they need to interact with, we simulate merging their coulomb potentials'''
def apply_merge_error(circuit, register, errors: dict):
    p = errors['merge'].p
    if p > 0:
        circuit.append(errors['merge'].op, register, p)

''' apply_split_error
Simulating splitting the coulomb potentials of check qubit and data qubit modules by appling an error probability to them.'''
def apply_split_error(circuit, register, errors: dict):
    p = errors['split'].p
    if p > 0:
        circuit.append(errors['split'].op, register, p)
