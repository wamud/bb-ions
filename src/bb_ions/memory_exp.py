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
        self.l = code.l
        self.m = code.m
        
        shape = (3 + (not reuse_check_qubits), self.l, self.m)
        
        arr_func = partial(self.gvw_to_k, self.l, self.m) 
        self.index_arr = np.fromfunction(arr_func, shape, dtype=int) 
        
        self.X = self.index_arr[0, :, :]
        self.L = self.index_arr[1, :, :]
        self.R = self.index_arr[2, :, :]
        
        z_group = 0 if reuse_check_qubits else 3
        self.Z = self.index_arr[z_group, :, :]

        

    @cached_property
    def data_registers(self):
        return join(self.L, self.R)

    @staticmethod
    def gvw_to_k(l: int, m: int, group: int, v: int, w: int):
        """
        convert group, row, coloumn index to integer index
        """
        return  group * l * m + v * m + w


class MemoryExperiment:
    """
    Handles generating stim circuit for a memory experiment.
    """
    def __init__(
        self, code: BBCode, arch: Architecture, reuse_check_qubits=True,
        memory_basis = "Z", exclude_opposite_basis_detectors=True,
        sequential=True, extraction_cycles = None 
    ):
        """
        args
        ----
        """
        if memory_basis not in ("X", "Z"):
            raise ValueError("memory_basis not X or Z")
        else:
            self.memory_basis = memory_basis
        
        self.code = code
        self.arch = arch
        self.reg = Registers(code, reuse_check_qubits)

        self.Junion = sorted(set(code.left_pow[1] + code.right_pow[1]))
        self.jval = self.Junion[0]
        self.exclude_opposite_basis_detectors = exclude_opposite_basis_detectors
        self.sequential = sequential
        self.extraction_cycles = extraction_cycles

    def idle_data_registers(self, idle_type, scaling=1, reg=None):
        
        if reg is None:
            targ = self.reg.data_registers
        else:
            targ = reg
            
        return self.arch.apply_probabilistic_gate(
            idle_type,
            targ,
            scaling=scaling
        )

    def apply_aligned_gates(
        self, gate: str, poly_pow: Iterable[Iterable[int]],
        control: np.ndarray, target: np.ndarray
    ):
        """
        Apply gates between control check_registers and target left / right
        data registers.
        """
        l = self.code.l
        m = self.code.m
        
        other_data_reg = self.reg.L if target is self.reg.R else self.reg.R
        circ = stim.Circuit()
        for i, j in zip(*poly_pow):
            if j == self.jval:
                #circ.append("S", self.jval + max(self.Junion) )
                #circ.append("Y", (i + max(self.Junion), j + max(self.Junion)))
                #circ.append("Z", target.flatten())
                for v in range(l):
                    for w in range(m):
                        circ += self.arch.apply_gate(
                            gate,
                            (control[v, w], target[(v + i) % l, (w + j) % m])
                       )
                    if self.sequential:
                        circ += self.idle_data_registers(
                            "idle_2q", reg=other_data_reg
                        )
                        for w in range(m):
                            for vprime in range(l):

                                if vprime != (v + i) % l:
                                    circ += self.idle_data_registers(
                                        "idle_2q", reg=target[vprime, w]
                                    )
                        circ.append("TICK")
                
                if not self.sequential:
                    circ += self.idle_data_registers("idle_2q", reg=other_data_reg)
                    circ.append("TICK")
        return circ

    def apply_cyclic_shift(self, basis):
        """
        """
        if basis not in ("X", "Z"):
            raise ValueError("basis not X or Z")

        check_register = getattr(self.reg, basis)
        if basis == "X":
            union_set = self.Junion
        else:
            union_set = tuple(reversed([-i for i in self.Junion]))

        #print("union_set ", union_set)
        #print(self.code.left_pow, self.code.right_pow)
        if basis == "X":
            gate = "CX"
            left_pow = self.code.left_pow
            right_pow = self.code.right_pow
        else:
            gate = "CZ"
            right_pow = [[-i for i in j] for j in self.code.left_pow]
            left_pow = [[-i for i in j] for j in self.code.right_pow]
        #print("left_pow", left_pow)
        #print("right_pow", right_pow)

        circ = stim.Circuit()
        #circ.append("X", [abs(a) for a in union_set])
        for j_idx, jval in enumerate(union_set):
            scaling = abs((jval % self.code.m) - (self.jval % self.code.m))

            circ += self.arch.apply_probabilistic_gate(
                    "shift", check_register, scaling=scaling
            )

            # Update to next configuration (last one looks uneeded)
            self.jval = jval
            
            #print(basis)
            #print(self.arch.apply_probabilistic_gate(
            #       "shift", check_register, scaling=scaling
            #))

            circ += self.idle_data_registers("idle_shift", scaling=scaling)
            #print(self.idle_data_registers("idle_shift", scaling=scaling))
            

            if scaling > 0:
                circ.append("TICK")
            
            circ += self.arch.apply_probabilistic_gate("shuttle", check_register)
            #print(self.arch.apply_probabilistic_gate("shuttle", check_register))
            circ += self.idle_data_registers("idle_shuttle")
            #print(self.idle_data_registers("idle_shuttle"))
            circ.append("TICK")
            
            circ += self.arch.apply_probabilistic_gate(
                "merge",
                join(check_register, self.reg.L, self.reg.R)
            )
            circ.append("TICK")
            

            #circ.append("X", 1)
            circ += self.apply_aligned_gates(
                gate, left_pow, check_register, self.reg.L
            )
            #print(basis)
            #print("\nLeft\n")
            #print(self.apply_aligned_gates(
            #    gate, left_pow, check_register, self.reg.L
            #))
        
            #circ.append("X", 2)
            circ += self.apply_aligned_gates(
                gate, right_pow, check_register, self.reg.R
            )
            #print("\n right \n")
            #print(self.apply_aligned_gates(
            #    gate, right_pow, check_register, self.reg.R
            #))

            ##circ.append("X", 3)

            circ += self.arch.apply_probabilistic_gate(
                "split",
                join(check_register, self.reg.L, self.reg.R)
            )
            circ.append("TICK")

            circ += self.arch.apply_probabilistic_gate("shuttle", check_register)
            circ += self.idle_data_registers("idle_shuttle")
            circ.append("TICK")
            

        return circ

    def apply_stabilizer_half_round(self, basis, round_0=False):
        """
        Perform stabilizer checks according to basis (either X or Z)
        """
        if basis not in ("X", "Z"):
            raise ValueError("basis not X or Z")

        check_register = getattr(self.reg, basis)
        
        # Initialise check register
        circ = self.arch.apply_gate("RZ", check_register)
        #circ.append(basis, (not round_0) and basis == "Z")
        
        # Idle data registers if not the first round X check
        if not round_0 or basis == "Z":
            circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")
        
        # Hadamard to |+> states
        # why not RX initialise?
        circ += self.arch.apply_gate("H", check_register)
        # Idle data registers if not the first round X check
        if not round_0 or basis == "Z":
            circ += self.idle_data_registers("idle_1q")
       
        # If round 0, initialise data registers into memory basis.
        # This is done before the X checks for the mesurement round
        if round_0 and basis == "X":
            circ += self.arch.apply_gate(
                f"R{self.memory_basis}", self.reg.data_registers
            )
        else:
            self.idle_data_registers("idle_1q")

        circ.append("TICK")
        
        circ += self.apply_cyclic_shift(basis)
        
        # Hadamard
        circ += self.arch.apply_gate("H", check_register)
        circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")

        # Measure
        circ += self.arch.apply_measurement("Z", check_register)
        circ += self.idle_data_registers("idle_meas")
        circ.append("TICK")
        
        n = self.code.qubits.physical
        if round_0 and self.memory_basis == basis:
            for i in range(n // 2, 0, -1):
                circ.append("DETECTOR", [stim.target_rec(-i)])
        elif not round_0:
            if (self.memory_basis == basis or 
                not self.exclude_opposite_basis_detectors):
                for i in range(n // 2, 0, -1):
                    circ.append(
                        "DETECTOR",
                        [stim.target_rec(-i), stim.target_rec(-i - n)]
                    )
        return circ

    def apply_stabilizer_round(self, round_0=False):
        """
        Apply a round of X and Z checks
        """
        circ = self.apply_stabilizer_half_round(basis="X", round_0=round_0)
        circ += self.apply_stabilizer_half_round(basis="Z", round_0=round_0)

        return circ

    def add_final_detectors(self):
        """
        add final detectors to X or Z checks depending on memory basis.
        """
        
        circ = stim.Circuit()
        check_mat = getattr(self.code.check_operators, self.memory_basis)
        checks = (np.nonzero(row)[0] for row in check_mat)
        
        n = self.code.qubits.physical
        for idx, c in enumerate(checks):
            if self.memory_basis == "Z":
                check_qubit_index  = -(n + n // 2) + idx 
            else:
                check_qubit_index = -2 * n + idx
        
            circ.append(
                "DETECTOR",
                [stim.target_rec(check_qubit_index)]
                + [stim.target_rec(q - n) for q in c]
            )

        return circ

    def add_logical_observables(self):
        """
        returns stim circuit to measure logical observable
        """
        circ = stim.Circuit()
        logical_op = getattr(self.code.logical_operators, self.memory_basis)
        n = self.code.qubits.physical
        record_indices = [np.nonzero(l)[0] - n for l in logical_op]
        
        for record in record_indices:
            circ.append(
                "OBSERVABLE_INCLUDE", [stim.target_rec(r) for r in record], 0.0
            )


        return circ

    @cached_property
    def stim_circuit(self):
        """
        Generate stim circuit for the memory experiment
        """

        circ = self.apply_stabilizer_round(round_0=True)
        repeated_loop = self.apply_stabilizer_round()

        if self.extraction_cycles is None:
            self.extraction_cycles = self.code.estimate_distance()
        circ += (self.extraction_cycles - 1) * repeated_loop

        # measure all data qubits
        circ += self.arch.apply_measurement(
            self.memory_basis, self.reg.data_registers
        )

        # Add final detectors
        circ += self.add_final_detectors()

        circ += self.add_logical_observables()

        circ.detecting_regions()
        
        return circ

