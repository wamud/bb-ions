from functools import partial
import numpy as np
from bb_ions.code import BBCode
from bb_ions.architecture import Architecture


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
        
        self.data_registers = np.concatenate((self.L, self.R), axis=None)

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
            memory_basis = "Z", exclude_opposite_basis_detectors=True):
        self.code = code
        self.arch = arch
        self.reg = Registers(code, reuse_check_qubits)
        
        self.Junion = sorted(set(self.left_pow[1] + self.right_pow[1]))
        self.jval = self.Junion[0]
        
        if memory_basis not in ("X", "Z"):
            raise ValueError("memory_basis not X or Z")
        else:
            self.memory_basis = memory_basis
    
    def initialise_registers(self):
        pass

    def idle_data_registers(self, idle_type, scaling=1):
        return self.arch.apply_probabilistic_gate(
                idle_type,
                self.reg.data_registers,
                scaling=scaling
        )

    def apply_cyclic_shift(self, check_register):
        """
        """
        circ = stim.Circuit()
        for jval in self.Junion:
            scaling = (jval % self.m) - (self.jval % self.m)
            circ += arch.apply_probabilistic_gate(
                    "shift", check_register, scaling=scaling
            )
            circ += self.idle_data_registers("idle_shift", scaling=scaling)
            

            circ.append("TICK")
            
            circ += self.arch.apply_probabilistic_gate("shuttle", check_register)
            circ += self.idle_data_registers("idle_shuttle")
            circ.append("TICK")
            
            circ += arch.apply_probabilistic_gate(
                    "merge",
                    np.concatenate(
                        (check_register, self.reg.L, self.reg.R), axis=None
                    )
            )

            return circ

    def apply_stabilizer_half_round(self, basis, round_0=False):
        """
        Peform stabilizer checks according to basis (either X or Z)
        """
        if basis not in ("X", "Z"):
            raise ValueError("basis not X or Z")

        check_register = getattr(self.reg, basis)
        
        # Initialise check register
        circ = self.arch.apply_gate("RZ", check_register)

        # Idle data registers if not the first round
        if not round_0:
            circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")
        
        # Hadamard to |+> states
        # why not RX initialise?
        circ += self.arch.apply_gate("H", check_register)
        
        if round_0 and basis == "X":
            # Initialise data qubits into memory basis
            circ += self.arch.apply_gate(
                f"R{self.memory_basis}", self.data_registers
            )
        else:
            self.idle_data_registers("idle_1q")

        circ.append("TICK")
        
        circ += self.apply_cyclic_shift(check_register)
        
        # Hadamard
        circ += self.arch.apply_gate("H", check_register)
        circ += self.idle_data_registers("idle_1q")
        circ.append("TICK")

        # Measure
        circ += self.arch.apply_probabilistic_gate("MZ", check_register)
        circ += self.idle_data_registers("idle_meas")

        n = self.code.physical
        if round_0 and self.memory_basis == basis:
            for i in range(n // 2, 0, -1):
                circ.append("DETECTOR", [stim.target_rec(-i)])
        elif not round_0:
            if self.memory_basis == basis or not exclude_opposite_basis_detectors:
                for i in range(n // 2, 0, -1):
                    circ.append(
                        "DETECTOR",
                        [stim.target_rec(-i), stim.target_rec(-i - n]
                    )
        
    def apply_stabilizer_round(self, round_0=False):
        """
        Apply a round of X and Z checks
        """

    def generate_stim_circuit(self):
        pass
