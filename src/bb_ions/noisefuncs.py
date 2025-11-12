''' noisefuncs'''

import numpy as np


LEG_SPACING = 5e-3 # [m]

SHUTTLE_DISTANCE = 2.5e-3 # 2.5mm to get in and out of a leg. Calling this a "shuttle" as opposed to the shuttling that occurs during a cyclic shift (which we call a "shift")

SHUTTLE_SPEED = 1 # [m/s]


class Error:
    def __init__(self, operation, p):
        self.op = operation # the error operation
        self.p = p # its probability



''' tham_modules_errors
Defines noise values as per Tham ... Delfosse "qubit modules" [2508.01879] (page 4) which uses Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module (i.e. assuming each module is a long chain) plus a cyclic shift error rate of 30p/100 when shifting modules (to align them).
(Note we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)'''
def tham_modules_errors(p):

    errors = {

        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 10), # note is depolarize (as per longcahin paper) as opposed to X_ERROR
        "RX" : Error("DEPOLARIZE1", p / 10), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", p / 10),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p / 10),
        "MX" : Error("Z_ERROR", p / 10),


        # "Qubit modules" [2508.01879] pg. 4
        "shift" : Error("DEPOLARIZE1", 30 * p / 100),
        "shift_const" : None, # shift_const is used to make errors proportional to the length of the shift. Tham et al. say the noise is independent of the length of the shift so set shift_const to None. (This means circfuncs.update_shift_probs will not change the shift error, making it independent of the length of the shift).


        # Additional for our architecture:
        "shuttle" : Error("DEPOLARIZE1", 0), # we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),

    }
    return errors

''' tham_modules_idle_errors
    Defines noise values as per Tham ... Delfosse "qubit modules" [2508.01879] (page 4) which uses Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module (i.e. assuming each module is a long chain) plus a cyclic shift error rate of 30p/100 when shifting modules (to align them).'''
def tham_modules_idle_errors(p):

    idle_during = {
        
        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 100),
        "RX" : Error("DEPOLARIZE1", p / 100), 
        "H" : Error("DEPOLARIZE1", p / 100),
        "CNOT" : Error("DEPOLARIZE1", p / 100),
        "CZ" : Error("DEPOLARIZE1", p / 100),
        "MZ" : Error("DEPOLARIZE1", 30 * p / 100),
        "MX" : Error("DEPOLARIZE1", 30 * p / 100),

        # "Qubit modules" [2508.01879] pg. 4
        "shift" : Error("DEPOLARIZE1", 30 * p / 100), 
        "shift_const" : None, # i.e. set to None means the the shift idling error will always be the value in the line above rather than this constant multiplied by the length of the shift

        # Additional for our architecture:
        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "pause" : Error("DEPOLARIZE1", 0), # this is an idling error applied at the beginning of each round of stabiliser measurements; simulates waiting before each round of stab. measurement

    }

    return idle_during


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

        "shift_const" : 0, 
        
    }
    
    return errors



def zero_idling():
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
        "pause" : Error("DEPOLARIZE1", 0),

        "shift_const" : 0,
    }

    return idle_during


''' uniform_errors
Defines a standard depolarising noise channel, as used in original Bravyi et al. BB paper [2308.0791] (see page 16).'''
def uniform_errors(p):

    errors = {

        "RZ" : Error("X_ERROR", p), 
        "RX" : Error("Z_ERROR", p),
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p),
        "MX" : Error("Z_ERROR", p),


        # Qubit module errors -- None
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_const" : None, 

        # Additional for our architecture - None
        "shuttle" : Error("DEPOLARIZE1", 0), 
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
    }
    
    return errors


def uniform_idling(p):

    idle_during = {
        "RZ" : Error("DEPOLARIZE1", p),
        "RX" : Error("DEPOLARIZE1", p), 
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE1", p),
        "CZ" : Error("DEPOLARIZE1", p),
        "MZ" : Error("DEPOLARIZE1", p),
        "MX" : Error("DEPOLARIZE1", p),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_const" : 0,

        "pause" : Error("DEPOLARIZE1", 0),

    }

    return idle_during


''' p_idle_dephasing
If a qubit is idling for time t, this function returns what the probability of a Z-error occurring on it will be. It takes as inputs t, the time the qubit is idling for, and T2, the characteristic time for dephasing of an idling qubit. For example p_idle(t = 100e-6, T2 = 50s) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience a Z-error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2))'''
def p_idle_dephasing(t, T2):
    p = 0.5 * (1 - np.exp(-t / T2))
    return p

''' dephasing_idle_errors
Defines idle errors as the realistic dephasing noise for a given T2 value with times set by the sources listed below in comments'''
def dephasing_idle_errors(T2):

    t_1q = 7.5e-6     # [27] in Bruzewicz et al. Table 1
    t_2q = 100e-6     # [27] in Bruzewicz et al. Table 1
    t_m  = 200e-6     # Myerson et al. https://doi.org/10.1103/PhysRevLett.100.200502
    t_r  = 207.5e-6   # reset is measurement + single qubit X
    
    idle_during = {
        
        "RZ" : Error("Z_ERROR", p_idle_dephasing(t_r, T2)),
        "RX" : Error("Z_ERROR", p_idle_dephasing(t_r, T2)), 
        "H" : Error("Z_ERROR", p_idle_dephasing(t_1q, T2)),    
        "CNOT" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)), 
        "CZ" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)),   
        "MZ" : Error("Z_ERROR", p_idle_dephasing(t_m, T2)),
        "MX" : Error("Z_ERROR", p_idle_dephasing(t_m, T2)),

        "shuttle" : Error("Z_ERROR", p_idle_dephasing(SHUTTLE_DISTANCE / SHUTTLE_SPEED, T2)), 
        
        "merge" : Error("Z_ERROR", 0),
        "split" : Error("Z_ERROR", 0),
        
        "shift" : Error("Z_ERROR", 0.1), # will be updated by circfuncs.update_shift_prob as long as p != 0.

        "shift_const" : T2, # T2 time

        "pause" : Error("Z_ERROR", 0),
    }

    return idle_during


def our_uniform_plus_shift_and_shuttle(p, T2):

    SHUTTLE_SPEED = 1 # [m/s]

    errors = {
        
        # Uniform errors:
        "RZ" : Error("X_ERROR", p), 
        "RX" : Error("Z_ERROR", p),
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p),
        "MX" : Error("Z_ERROR", p),

        # Qubit module errors
        "shift_const" : T2, 
        "shift" : Error("Z_ERROR", 0.1), # will be updated each shift
 

        # Additional for our architecture 
        "shuttle" : Error("Z_ERROR", p_idle_dephasing(SHUTTLE_DISTANCE / SHUTTLE_SPEED, T2)), 

        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
    }
    
    return errors






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
