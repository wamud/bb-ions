from functools import partial, cached_property
from typing import Iterable
import numpy as np
import stim
from bb_ions.code import BBCode
from bb_ions.architecture import Architecture


def join(*args):
    """
    Return the provided numpy arrays as a single 1D array.
    """
    return np.concatenate(args, axis=None)


class Registers:
    def __init__(self, code: BBCode, reuse_check_qubits: bool = True):
        """
        Class to manage register indices of X-check, Z-check and L (left) and
        R (right) data qubits. Registers can be accessed as 2D arrays of
        dimension (l, m).
        """
        # BB code parameters
        self.l = code.l
        self.m = code.m

        # Shape of total register array (increases dimension if different
        # registers are used for Z and X checks.
        shape = (3 + (not reuse_check_qubits), self.l, self.m)

        # Function to generate total register array
        arr_func = partial(self.gvw_to_k, self.l, self.m)

        # Generate total register array
        self.index_arr = np.fromfunction(arr_func, shape, dtype=int)

        # Slice X, L and R arrays
        self.X = self.index_arr[0, :, :]
        self.L = self.index_arr[1, :, :]
        self.R = self.index_arr[2, :, :]

        # Slice Z-check array depending on whether using same registers as X.
        z_group = 0 if reuse_check_qubits else 3
        self.Z = self.index_arr[z_group, :, :]

    @cached_property
    def data_registers(self) -> np.ndarray:
        """
        Return data registers (L + R) as a 1D array.
        """
        return join(self.L, self.R)

    @staticmethod
    def gvw_to_k(l: int, m: int, group: int, v: int, w: int) -> int:
        """
        convert group, row, coloumn index to integer index
        """
        return group * l * m + v * m + w


class MemoryExperiment:
    """
    Handles generating stim circuit for a memory experiment.
    """

    def __init__(
        self,
        code: BBCode,
        arch: Architecture,
        reuse_check_qubits=True,
        memory_basis="Z",
        exclude_opposite_basis_detectors=True,
        sequential=True,
        extraction_cycles=None,
    ):
        """
        args
        ----
        code: BBCode object for a given BB code.
        arch: Architecture object for a specific architecture.
        reuse_check_qubits: flags whether to use same registers for X and Z
            checks.
        memory_basis: "X" or "Z" determining the memory basis of the
            experiment.
        exclude_opposite_basis_detectors: flag to exclude detectors from basis
            other than memory_basis.
        sequential: Whether to apply 2-qubit gates sequentially within a module.
        extraction_cycles: Number of extraction cycles, defaults to max distance
            determined by bposd.
        """
        if memory_basis not in ("X", "Z"):
            raise ValueError("memory_basis not X or Z")
        else:
            self.memory_basis = memory_basis

        self.code = code
        self.arch = arch
        self.reg = Registers(code, reuse_check_qubits)

        # All y base powers for X checks (Z basis ones are negation of this)
        self.Junion = sorted(set(code.left_pow[1] + code.right_pow[1]))
        # holds current y power configuration of the architecture.
        self.jval = self.Junion[0]
        self.exclude_opposite_basis_detectors = exclude_opposite_basis_detectors
        self.sequential = sequential
        self.extraction_cycles = extraction_cycles

    def idle_data_registers(self, idle_type, scaling=1, reg=None) -> stim.Circuit:
        """
        Idle data register for idle_type (key in arch.errors) with with
        error rate scaled by `scaling`. By default acts on all data registers
        but can be overriden by passing a reg argument to act on.
        """
        if reg is None:
            targ = self.reg.data_registers
        else:
            targ = reg

        return self.arch.apply_probabilistic_gate(idle_type, targ, scaling=scaling)

    def apply_aligned_gates(
        self,
        gate: str,
        poly_pow: Iterable[Iterable[int]],
        control: np.ndarray,
        target: np.ndarray,
    ) -> stim.Circuit:
        """
        Apply gates between control check_registers and target left / right
        data registers.

        args
        ----
        gate: key of gate defined in arch.errors.
        poly_pow: polynomial powers for block (A, B, A^T, B^T)
        control: check registers (X or Z)
        target: data registers (L or R)
        """
        # Code parameters
        l = self.code.l
        m = self.code.m

        # Data register not being acted on L / R.
        other_data_reg = self.reg.L if target is self.reg.R else self.reg.R

        # Initialise stim circuit.
        circ = stim.Circuit()

        # Loop over all polynomial powers for this block and apply CX / CZ
        # gates for ancillas aligned with data modules.
        for i, j in zip(*poly_pow):
            # Check if current power is aligned.
            if j == self.jval:
                # Apply gate for each aligned qubit in module
                for v in range(l):
                    for w in range(m):
                        circ += self.arch.apply_gate(
                            gate, (control[v, w], target[(v + i) % l, (w + j) % m])
                        )
                    if self.sequential:
                        # Idle other data register after applying one gate
                        # within each module.
                        circ += self.idle_data_registers("idle_2q", reg=other_data_reg)

                        # Idle non-participating qubits in current data register.
                        for w in range(m):
                            for vprime in range(l):
                                if vprime != (v + i) % l:
                                    circ += self.idle_data_registers(
                                        "idle_2q", reg=target[vprime, w]
                                    )
                        circ.append("TICK")

                # Apply gates within a module parallelly
                if not self.sequential:
                    circ += self.idle_data_registers("idle_2q", reg=other_data_reg)
                    circ.append("TICK")
        return circ

    def apply_cyclic_shift(self, basis: str) -> stim.Circuit:
        """
        Apply a cyclic shifts and perform ancillae checks for a given basis.
        Will do CNOT checks for A and B if basis is "X" and CZ checks for B^T
        A^T if basis is "Z" cycling ancillae qubits appropriately.

        args
        ----
        basis: "X" or "Z".
        """
        if basis not in ("X", "Z"):
            raise ValueError("basis not X or Z")

        # Get X or Z registers by field name
        check_register = getattr(self.reg, basis)
        # Define y powers to cycle through
        if basis == "X":
            union_set = self.Junion
        else:
            # Reversing to keep compatibility with earlier code.
            union_set = tuple(reversed([-i for i in self.Junion]))

        # Define gates and polynomial powers to cycle through for the current
        # basis.
        if basis == "X":
            gate = "CX"
            left_pow = self.code.left_pow
            right_pow = self.code.right_pow
        else:
            gate = "CZ"
            right_pow = [[-i for i in j] for j in self.code.left_pow]
            left_pow = [[-i for i in j] for j in self.code.right_pow]

        # Initialise stim circuit.
        circ = stim.Circuit()

        # Loop over y powers and apply gates for the current block.
        for j_idx, jval in enumerate(union_set):
            # Calculate shift distance to scale shift error and idle time.
            scaling = abs((jval % self.code.m) - (self.jval % self.code.m))

            # Apply shift error
            circ += self.arch.apply_probabilistic_gate(
                "shift", check_register, scaling=scaling
            )

            # Update current y power configuration
            self.jval = jval

            # Apply shift idling error on data registers
            circ += self.idle_data_registers("idle_shift", scaling=scaling)

            # If scaling == 0, no shift error is applied so time step doesn't
            # exist.
            if scaling > 0:
                circ.append("TICK")

            # Apply shuttling error (racetrack -> leg)
            circ += self.arch.apply_probabilistic_gate("shuttle", check_register)

            # Idle data registers during shuttling.
            circ += self.idle_data_registers("idle_shuttle")
            circ.append("TICK")

            # Apply merge error
            circ += self.arch.apply_probabilistic_gate(
                "merge", join(check_register, self.reg.L, self.reg.R)
            )
            circ.append("TICK")

            # Apply aligned check gates between check register and left data
            # qubits
            circ += self.apply_aligned_gates(gate, left_pow, check_register, self.reg.L)

            # Apply aligned check gates between check register and right data
            # qubits
            circ += self.apply_aligned_gates(
                gate, right_pow, check_register, self.reg.R
            )

            # Apply split error
            circ += self.arch.apply_probabilistic_gate(
                "split", join(check_register, self.reg.L, self.reg.R)
            )
            circ.append("TICK")

            # Apply shuttling error (leg -> racetrack)
            circ += self.arch.apply_probabilistic_gate("shuttle", check_register)

            # Idle data registers during shuttling
            circ += self.idle_data_registers("idle_shuttle")
            circ.append("TICK")

        return circ

    def apply_stabilizer_half_round(self, basis, round_0=False) -> stim.Circuit:
        """
        Perform stabilizer checks according to basis (either X or Z)

        args
        ----
        basis: "X" or "Z"
        round_0: flags whether it is the first round of stabilizer checks which
            involve separate operations such as data register initialisation
        """
        if basis not in ("X", "Z"):
            raise ValueError("basis not X or Z")

        # Get X or Z registers by field name
        check_register = getattr(self.reg, basis)

        # Initialise check register
        circ = self.arch.apply_gate("RZ", check_register)

        # Idle data registers if not the first round X check
        if not round_0 or basis == "Z":
            circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")

        # Hadamard to |+> states
        circ += self.arch.apply_gate("H", check_register)

        # Idle data registers if not the first round X check
        if not round_0 or basis == "Z":
            circ += self.idle_data_registers("idle_1q")

        # If round 0, initialise data registers into memory basis.
        # This is done before the X checks for the measurement round
        if round_0 and basis == "X":
            circ += self.arch.apply_gate(
                f"R{self.memory_basis}", self.reg.data_registers
            )
        else:
            self.idle_data_registers("idle_1q")

        circ.append("TICK")

        # Apply cyclic shifts and perform CX / CZ gates for the current basis.
        circ += self.apply_cyclic_shift(basis)

        # Hadamard check registers
        circ += self.arch.apply_gate("H", check_register)

        # Idle data registers during Hadamard
        circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")

        # Measure check registers
        circ += self.arch.apply_measurement("Z", check_register)

        # Idle data registers during measurement.
        circ += self.idle_data_registers("idle_meas")
        circ.append("TICK")

        n = self.code.qubits.physical
        if round_0 and self.memory_basis == basis:
            for i in range(n // 2, 0, -1):
                # add detectors to memory basis checks in round 0.
                circ.append("DETECTOR", [stim.target_rec(-i)])
        elif not round_0:
            if self.memory_basis == basis or not self.exclude_opposite_basis_detectors:
                for i in range(n // 2, 0, -1):
                    # add detectors to memory basis checks or all checks
                    # depending on exclude_opposite_basis_detectors.
                    circ.append(
                        "DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - n)]
                    )
        return circ

    def apply_stabilizer_round(self, round_0=False) -> stim.Circuit:
        """
        Apply a round of X and Z checks

        args
        ----
        round_0: flag first round
        """
        # Apply X round
        circ = self.apply_stabilizer_half_round(basis="X", round_0=round_0)

        # Apply Z round
        circ += self.apply_stabilizer_half_round(basis="Z", round_0=round_0)

        return circ

    def add_final_detectors(self) -> stim.Circuit:
        """
        add final detectors to X or Z checks depending on memory basis.
        """

        # Initialise stim circuit.
        circ = stim.Circuit()

        # Get check matrix based on `self.memory_basis`.
        check_mat = getattr(self.code.check_operators, self.memory_basis)

        # Get non-zero entries for each row of the check matrix.
        checks = (np.nonzero(row)[0] for row in check_mat)

        n = self.code.qubits.physical

        # Add detectors based on memory basis
        for idx, c in enumerate(checks):
            if self.memory_basis == "Z":
                check_qubit_index = -(n + n // 2) + idx
            else:
                check_qubit_index = -2 * n + idx

            circ.append(
                "DETECTOR",
                [stim.target_rec(check_qubit_index)]
                + [stim.target_rec(q - n) for q in c],
            )

        return circ

    def add_logical_observables(self):
        """
        returns stim circuit to measure logical observable
        """
        # Initialise stim circuit.
        circ = stim.Circuit()

        # Get logical op based on self.memory_basis.
        logical_op = getattr(self.code.logical_operators, self.memory_basis)
        n = self.code.qubits.physical

        # Get non-zero entries in logical operator per row
        record_indices = [np.nonzero(lg)[0] - n for lg in logical_op]

        # Add detectors for logical observable
        for record in record_indices:
            circ.append("OBSERVABLE_INCLUDE", [stim.target_rec(r) for r in record], 0.0)

        return circ

    @cached_property
    def stim_circuit(self):
        """
        Generate stim circuit for the memory experiment. This is a cached
        property which is computed on first access and cached value returned
        thereafter.
        """

        # Generate round 0 stim circuit
        circ = self.apply_stabilizer_round(round_0=True)

        # Generate repeated syndrome extraction circuit
        repeated_loop = self.apply_stabilizer_round()

        # Determine number of extraction cycles. Defaults to max distance
        # given by bposd.
        if self.extraction_cycles is None:
            self.extraction_cycles = self.code.estimate_distance()

        # Generate remaining syndrome rounds
        circ += (self.extraction_cycles - 1) * repeated_loop

        # measure all data qubits
        circ += self.arch.apply_measurement(self.memory_basis, self.reg.data_registers)

        # Add final detectors
        circ += self.add_final_detectors()

        # Add logical operators
        circ += self.add_logical_observables()

        # Test that detectors are valid.
        circ.detecting_regions()

        return circ
