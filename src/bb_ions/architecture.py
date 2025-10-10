from typing import Iterable
import stim


def default_errors(p=0.1):
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
        "idle" : ("DEPOLARIZE1", p / 100),
        "shuttle" : ("DEPOLARIZE1", p / 10),
        "merge" : ("DEPOLARIZE1", p / 10),
        "split" : ("DEPOLARIZE1", p / 10),
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
 
        gate_err, err_rate = self.errors[gate_name]
        
        circ = stim.Circuit()
        circ.append(gate_name, register)

        if err_rate > 0:
            circ.append(gate_err, register, err_rate)

        return circ

    def apply_probabilistic_gate(self, op_name: str, registers: Iterable[int]):
        """
        Apply operations with probabilistic effects such as idling, shuttling
        and measurement.
        """
        op_err, err_rate = self.errors[op_name]

        circ = stim.Circuit()
        circ.append(op_err, registers, err_rate)

        return circ

    def apply_shift(self, shift):
        pass
