""" 
Hardware-inspired noise and error functions relevant to our architecture
"""


class NoiseTimes:
    def __init__(
        self,
        t_init,
        t_had,
        t_merge,
        t_split,
        t_cnot,
        t_cz,
        t_shuttle,
        t_shift_const,
        t_meas,
        t_idle,
        t_idle_meas,
    ):
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


import numpy as np



def make_uniform_noisetimes(t):
    """ 
    Create an object noisetimes of class NoiseTimes to store the lengths of
    time each operation takes so they can be used when applying noise. This
    sets them all uniformly to the same input to this function t. The object
    noisetimes can be modified afterwards to set specific times
    """

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

    noisetimes = NoiseTimes(
        t_init,
        t_had,
        t_merge,
        t_split,
        t_cnot,
        t_cz,
        t_shuttle,
        t_shift_const,
        t_meas,
        t_idle,
        t_idle_meas,
    )

    return noisetimes



def make_longchain_noisetimes(t):
    """ 
    Create an object noisetimes of class NoiseTimes to store the lengths of
    time each operation takes so they can be used when applying noise. Sets
    according to Ye Delfosse longchain: 2503.2207
    """

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

    noisetimes = NoiseTimes(
        t_init,
        t_had,
        t_merge,
        t_split,
        t_cnot,
        t_cz,
        t_shuttle,
        t_shift_const,
        t_meas,
        t_idle,
        t_idle_meas,
    )

    return noisetimes



def make_uniform_agnostic_noisetimes(t):
    """ 
    Sets noise of all operations to the input t except removes the noise
    introduced by our hardware proposal, namely merge, split, shuttle, shift.
    """

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

    noisetimes = NoiseTimes(
        t_init,
        t_had,
        t_merge,
        t_split,
        t_cnot,
        t_cz,
        t_shuttle,
        t_shift_const,
        t_meas,
        t_idle,
        t_idle_meas,
    )

    return noisetimes



def make_longchain_agnostic_noisetimes(t):
    """ 
    Sets noise of all operations as per Ye Delfosse longchain paper 2503.2207.
    Also removes noise introduced by our hardware proposal, namely merge,
    split, shuttle, shift are set to zero.
    """

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

    noisetimes = NoiseTimes(
        t_init,
        t_had,
        t_merge,
        t_split,
        t_cnot,
        t_cz,
        t_shuttle,
        t_shift_const,
        t_meas,
        t_idle,
        t_idle_meas,
    )

    return noisetimes



def p_init(t):
    """ 
    The probability of initalising to an orthogonal eigenstate. Depends on time
    taken to initalise (t_init). todo: fill with more than dummy function.
    """

    p = t
    return p



def p_had(t):
    """
    The probability of an error occuring on the hadamard gate. For what error
    is actually applied see 'hadamard' in circfuncs.
    """

    p = t
    return p



def p_shuttle(t_shuttle):
    """ 
    If shuttling a qubit it will experience an error associated with
    acceleration, time of shuttle, number of T junctions it passes through and
    deceleration. todo: make this function more than a dummy
    """

    p = t_shuttle
    return p



def p_shift(t_shift):
    """ 
    todo: make more than dummy function
    """

    p = t_shift
    return p



def p_merge(t):
    """ 
    Depends on length of time merging the coulomb potential of two ion modules
    takes. Todo: fill this function. At the moment is just a dummy function
    returning t.
    """

    p = t
    return p



def p_cnot(t_cnot):
    """
    todo: fill 
    """

    p = t_cnot
    return p



def p_cz(t_cz):
    """
    todo: fill 
    """

    p = t_cz
    return p



def p_split(t):
    """ 
    todo: make more than dummy function
    """

    p = t
    return p



def p_meas(t_meas):
    """ 
    Probability of a measurement error given the time of measurement was
    t_meas. At the moment just returns p equal to t_meas
    """

    p = t_meas
    return p




def p_idle(t_idle):
    """
    p_idle
    """
    p = t_idle
    return p


# ''' p_idle
# If a qubit is idling for time t, this function returns what the probability
# of an error occurring on it will be. This probability can be fed into other
# functions to say what the error will actually be. For example p_idle(T =
# 100e-6) = 1e-6 , indicating that if a qubit is idling for 100μs it will
# experience an error with probability 1e-6. This function assumes idling is a
# dephasing noise channel with p = 0.5(1 - e^(-t/T_2)) where T_2 = 50
# seconds'''
# def p_idle(t):
#     T2 = 50
#     p = 0.5 * (1 - np.exp(-t / T2))
#     return p



def apply_shuttle_error(circuit, register, t_shuttle):
    """ 
    Applies a depolarising noise channel of strength p to qubits in 'register'
    and in stim circuit 'circuit'. todo: make more accurate noise model once we
    have the info
    """

    p = p_shuttle(t_shuttle)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)


# ''' idle
# Adds dephasing (Z) noise of strength p (a Z operation is applied with probability p) to qubits in 'register.'''
# def idle(circuit, register, t_idle = 0):
#   p = p_idle(t_idle)
#   if p > 0:
#     circuit.append("Z_ERROR", register, p)



def idle(circuit, register, t_idle=0):
    """ 
    Adds depolarising noise of strength p_idle(t_idle) to qubits in register
    """

    p = p_idle(t_idle)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)



def apply_shift_error(circuit, register, t_shift):
    """ 
    todo: make more than dummy function
    """

    p = p_shift(t_shift)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)



def apply_merge_error(circuit, register, t_merge):
    """ 
    Once check modules have been cyclically shifted to the data qubit module
    they need to interact with, we simulate merging their coulomb potentials
    """

    p = p_merge(t_merge)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)



def apply_split_error(circuit, register, t_split):
    """ 
    Simulating splitting the coulomb potentials of check qubit and data qubit
    modules by appling an error probability to them. At the moment is just
    depolarising noise of strength p = t_split
    """

    p = p_split(t_split)
    if p > 0:
        circuit.append("DEPOLARIZE1", register, p)
