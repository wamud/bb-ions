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
    def __init__(self, code: BBCode, arch: Architecture, reuse_check_qubits=True):
        self.code = code
        self.reg = Registers(code, reuse_check_qubits)
    
    def initialise_registers(self):
        pass


    def generate_stim_circuit(self):
        pass

