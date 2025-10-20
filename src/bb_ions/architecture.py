from typing import Iterable
import stim
import numpy as np


def default_errors(p=0.01):
    """
    Defines default gates and error rates.
    """
    errors = {
        "RZ": ("DEPOLARIZE1", p / 10),
        "RX": ("DEPOLARIZE1", p / 10),
        "H": ("DEPOLARIZE1", p / 10),
        "CX": ("DEPOLARIZE2", p),
        "CZ": ("DEPOLARIZE2", p),
        "MZ": ("DEPOLARIZE1", p / 10),
        "MX": ("DEPOLARIZE1", p / 10),
        "idle_1q": ("DEPOLARIZE1", p / 100),
        "idle_2q": ("DEPOLARIZE1", p / 100),
        "idle_shift": ("DEPOLARIZE1", p / 100),
        "idle_shuttle": ("DEPOLARIZE1", p / 100),
        "idle_meas": ("DEPOLARIZE1", 30 * p / 100),
        "shuttle": ("DEPOLARIZE1", p / 10),
        "merge": ("DEPOLARIZE1", p / 10),
        "split": ("DEPOLARIZE1", p / 10),
        "shift": ("DEPOLARIZE1", p / 10),
    }

    return errors


class Architecture:
    """
    Class to handle architecture specific stim subcircuits
    """

    def __init__(self, errors=None):
        """
        args
        ----
        errors: dictionary holding key value pairs
            {op_name: (stim_gate_name / error, error_rate)}
        """
        if errors is None:
            errors = default_errors()
        self.errors = errors

    def apply_gate(
        self, gate_name: str, register: Iterable[int] | np.ndarray
    ) -> stim.Circuit:
        """
        Returns the stim circuit applying the gate and correspodning error
        in self.errors to the provided registers.

        args
        ----
        gate_name: key value in self.errors
        register: registers to apply gate to
        """

        # Flatten 2d arrays if needed.
        if isinstance(register, np.ndarray) and register.ndim > 1:
            register = register.flatten()

        # Get gate error type and error rate
        gate_err, err_rate = self.errors[gate_name]

        # Initalise a stim circuit
        circ = stim.Circuit()

        # Apply stim gate to registers
        circ.append(gate_name, register)

        # Apply gate error process if error rate is non-zero
        if err_rate > 0:
            circ.append(gate_err, register, err_rate)

        return circ

    def apply_probabilistic_gate(
        self, op_name: str, register: Iterable[int], scaling: int = 1
    ) -> stim.Circuit:
        """
        Apply operations with probabilistic error effects such as idling,
        shuttling. Allows error rate to be scaled by `scaling`.

        args
        ----
        op_name: key value in self.errors
        register: registers to apply operation to
        scaling: scaling factor for the error rate
        """

        # Get stim operation and error rate from self.errors
        op_err, err_rate = self.errors[op_name]

        # scale error rate (useful for different shift sizes)
        err_rate *= scaling

        # Flatten 2d arrays if needed.
        if isinstance(register, np.ndarray) and register.ndim > 1:
            register = register.flatten()

        # Initialise stim circuit
        circ = stim.Circuit()

        # Apply operation if error rate is non-zero.
        if err_rate > 0:
            circ.append(op_err, register, err_rate)

        return circ

    def apply_measurement(
        self, basis: str, register: Iterable[int] | np.ndarray
    ) -> stim.Circuit:
        """
        return stim circuit to measure register in given basis

        args
        ----
        basis: "X" or "Z".
        """

        meas_str = f"M{basis}"
        meas_err, err_rate = self.errors[meas_str]

        # Flatten 2d arrays if needed.
        if isinstance(register, np.ndarray) and register.ndim > 1:
            register = register.flatten()

        circ = stim.Circuit()
        # apply error before measurement (as opposed to gate)
        if err_rate > 0:
            circ.append(meas_err, register, err_rate)

        # apply measurement
        circ.append(meas_str, register)

        return circ
