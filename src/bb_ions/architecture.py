from typing import Iterable
import math
import stim
import numpy as np

def default_errors(p=0.01):
    """
    Defines default gates and error rates.
    """
    errors = {
        "RZ" : ("DEPOLARIZE1", p / 10),
        "RX" : ("DEPOLARIZE1", p / 10),
        "H" : ("DEPOLARIZE1", p / 10),
        "CX" : ("DEPOLARIZE2", p),
        "CZ" : ("DEPOLARIZE2", p),
        "MZ" : ("MZ", p / 10),
        "MX" : ("MX", p / 10),
        "idle_1q" : ("DEPOLARIZE1", p / 100),
        "idle_2q" : ("DEPOLARIZE1", p / 100),
        "idle_meas" : ("DEPOLARIZE1", 30 * p / 100),
        "shuttle" : ("DEPOLARIZE1", p / 10),
        "merge" : ("DEPOLARIZE1", p / 10),
        "split" : ("DEPOLARIZE1", p / 10),
        "shift" : ("DEPOLARIZE1", p / 10),
        
    }

    return errors


class Architecture:
    """
    Class to handle architecture specific stim subcircuits
    """
    def __init__(self, errors: dict = default_errors(), gate_times: dict = {}):
        self.errors = errors
        self.gate_times = gate_times
    
    def apply_gate(self, gate_name: str, register: Iterable[int]):
        """
        Returns the stim circuit applying the gate and error
        """
        
        # Flatten 2d arrays if needed.
        if isinstance(register, np.ndarray) and register.ndim > 1:
            register = register.flatten()

        gate_err, err_rate = self.errors[gate_name]
        
        circ = stim.Circuit()
        circ.append(gate_name, register)

        if err_rate > 0:
            circ.append(gate_err, register, err_rate)

        return circ

    def apply_probabilistic_gate(
            self, op_name: str, register: Iterable[int], scaling: int = 1):
        """
        Apply operations with probabilistic effects such as idling, shuttling
        and measurement.
        """
        op_err, err_rate = self.errors[op_name]
        # scale error rate (useful for different shift sizes)
        err_rate *= scaling

        # Flatten 2d arrays if needed.
        if isinstance(register, np.ndarray) and register.ndim > 1:
            register = register.flatten()
        
        circ = stim.Circuit()
        circ.append(op_err, register, err_rate)

        return circ
