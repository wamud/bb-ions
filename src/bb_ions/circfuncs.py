''' circfuncs
Given the parity check matrices (constructed using bbfuncs.py) and paramaters & logical operators (found using bbparams.py) of a Bicycle Bivariate [2308.07915] code, these functions (culminating in make_circuit at the bottom of this file) are for constructing a stim circuit that realises, with appropriate detectors annotated, a memory experiment using the BB code, i.e. it prepares logical |0⟩ or |+⟩ in all the logical qubits of the BB code, runs multiple rounds of stabiliser measurements, then measures all the data qubits.
Note Stim api reference: https://github.com/quantumlib/Stim/wiki/Stim-v1.9-Python-API-Reference'''

import stim
from .kfuncs import *
from .noisefuncs import *


class Registers:
    def __init__(self, C = None, L = None, R = None, X = None, Z = None, qL = None, qR = None, qC = None, qX = None, qZ = None):
        self.L = L
        self.R = R
        self.X = X
        self.Z = Z
        self.qL = qL
        self.qR = qR
        self.qX = qX
        self.qZ = qZ



''' make_registers
Makes lists of qubit indices dividing n = 2lm qubits evenly into qA, qB, qC, qD.
  qA: qubits 0 to n/2 - 1
  qB: qubits n/2 to n - 1
  qC: qubits n to 3n/2 - 1 
  qD: qubits 3n/2 to 2n - 1
This is the same as the make_registers function, it's just that it uses the convtok within it'''
def make_registers(l, m, reuse_check_qubits = True):
  
  
  if reuse_check_qubits == True:  # qX and qZ will be the same register
    X = 0
    L = 1
    R = 2
    Z = 0

  else:
    X = 0
    L = 1
    R = 2
    Z = 3

  qX = [convtok(l, m, X, v, w) for v in range(l) for w in range(m)]
  qL = [convtok(l, m, L, v, w) for v in range(l) for w in range(m)]
  qR = [convtok(l, m, R, v, w) for v in range(l) for w in range(m)]
  qZ = [convtok(l, m, Z, v, w) for v in range(l) for w in range(m)]
  
  registers = Registers(X = X, L = L, R = R, Z = Z, qX = qX, qL = qL, qR = qR, qZ = qZ)

  return registers

''' init
Sets qubits in the list 'register' to
- |0⟩ if basis == 'Z'
- |+⟩ if basis == 'X'
Also adds a depolarizing error (in alignment with lonchain paper 2503.2207) with probability p if 'longchain' in the noise parameter, otherwise the more usual reset error which sets to orthog. eigenstate'''
def init(basis, circuit, register, errors: dict):
  
  reset = f"R{basis}"
  circuit.append(reset, register) # append RZ or RX
  p = errors[reset].p

  if p > 0:
    error_op = errors[reset].op
    circuit.append(error_op, register, p)



''' hadamard
Appends a hadamard gate to the stim circuit on qubit(s) specified in 'register'.
After the gate it adds a depolarising noise of strength p (i.e. an error will occur with prob p. Given it occurs, pick one of X, Y or Z at random)'''
def hadamard(circuit, register, errors: dict):

  circuit.append("H", register)
  
  p = errors['H'].p
  
  if p > 0:
    circuit.append(errors['H'].op, register, p)


''' measure
Appends a 'basis'-basis measurement onto the qubits specified by 'register' on an input stim circuit 'circuit'. Probability of measurement error given by p_meas(t_meas)'''
def measure(basis: str, circuit, register, errors: dict):
  
  measure_string = f"M{basis}"
  
  p = errors[measure_string].p
  
  # Append error before measurement
  if p > 0: 
    error_op = errors[measure_string].op
    circuit.append(error_op, register, p)

  # Append measurement:
  circuit.append(measure_string, register)


''' tick
Appends a 'TICK' annotation to an input stim circuit, indicating the end of a time-step. '''
def tick(circuit):
  circuit.append("TICK")



''' add_final_detectors
After measuring all the data qubits in the X-basis (memory X) or Z-basis (memory Z) we want to check that each check qubit
has correctly reported the parity of the data qubits it was supposed to have measured.
Consequently, we add detectors that include each check qubit's parity multiplied with the parity of all the data qubits it checked.
Note that in mem. Z we only add these detectors to the Z-checks because the Z-measurements don't commute with X-stabiliser measurements (so we cannot check the parity of X-stabilisers from the data qubits that have been measured in Z
(whatever the parity may have been for an X-stab measurement, e.g. -1 for XXXX if |0000⟩ - |1111⟩ , once we collapse to, say, |1111⟩ there is no more X-parity check)).
Conversely for memory X we only add these detectors to the X-checks.

Optional further reading on inner workings of the function: figuring out each stim.target_rec[]:

In the final round, we performed X-check measurements, then Z-check measurements and finally all the data qubit measurements.
So the most recent rec (rec[-1]) is the n-th data qubit. It was the last measured. This is the last data qubit in qR. The zeroth data qubit, also the
zeroth data qubit in qL, is rec[-n]. Before this is qX then qZ check qubits.

So rec[ ]:
-1 to -n is data qubit measurements
-(n + 1) to - (n + n//2) is qZ measurements
-(n + n//2 + 1) to -2n is qX measurements

For each check qubit we want the detector to include its parity multiplied with the parity of all the data qubits it checked. These are contained in Hz and Hz, except the j-th data qubit will be rec[j - n] (so the 0-th data qubit was measured n measurements ago, the last / (n - 1)th data qubit was measured one measurement ago etc.'''
def add_final_detectors(circ, code, memory_basis):

  Hx = code.Hx
  Hz = code.Hz
  n = code.n

  if memory_basis == 'Z':
      
      for k in list(range(n//2)): # loop over Z-check qubits
          
          this_check_qubits_data_qubits = np.nonzero(Hz[k])[0] # extract the 1 positions (data qubit indices) in this check qubit's row of Hz

          circ.append("DETECTOR", 
          [stim.target_rec(-(n + n//2) + k)] # This Z-check qubit
          +
          [stim.target_rec(j - n) for j in this_check_qubits_data_qubits] # It's data qubits
          # Note it is j - n because if j is the last data qubit, j = n - 1, we need stim.target_rec(-1) (most recent measurement)
          )

  elif memory_basis == 'X':
      
      for k in list(range(n//2)):

          this_check_qubits_data_qubits = np.nonzero(Hx[k])[0]

          circ.append("DETECTOR", 
          [stim.target_rec(-2 * n + k)] # This X-check qubit
          +
          [stim.target_rec(j - n) for j in this_check_qubits_data_qubits] # Its data qubits
          )
  
  else:
    raise ValueError("Paramater 'memory_basis' must be either 'X' or 'Z' ")


'''get_nonzero_indices
For an array, this function returns a list (per row of the initial array) containing the indices of the nonzero terms.
'''
def get_nonzero_indices(array):

  array_indices = []

  for i in range(array.shape[0]):
    array_indices.append(np.nonzero(array[i])[0])

  return array_indices



''' add_logical_observables
The circuit ends with a parity check of the logical operators / observables. In a surface code, for example, XL and ZL are just vertical or horizontal chains of X's and Z's across the lattice.
In a BB code they are also chains of X's and Z's but on specific qubits, contained in the arrays Lx and Lz. For each logical qubit Lx and Lz contain a pair of anti-commuting logical operators.
These commute with the logical operators of other logical qubits. This function adds the Lx operators as observable if we are in memory X (preserving an eigenstate of Lx's, i.e. |+⟩_L)
or the Lz operators if we are in memory Z.

Inner workings (optional reading):
The data qubit measurements were the most recent measurements in the circuit so are from rec[-1] (the n-th data qubit) to rec[-n] (the zeroth data qubit). So if Lx has a 1 on qubit 0 for example, this is rec[-n]. Of a 1 on qubit 3 this means it needs rec[3 - n]. Just minus n from the index in Lx.'''
def add_logical_observables(circuit, n, Lx, Lz, memory):

  L = Lx if memory == 'X' else Lz

  num_logical_ops = L.shape[0]

  indices = get_nonzero_indices(L) # instead of L being 1's and 0's (like a parity check matrix) just make it a list of the indices of the 1's

  for i in range(num_logical_ops): # for each logical qubit

    recordings = (indices[i] - n).astype(int) # the measurements -- 'inner workings' note above

    circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(r) for r in recordings], 0.0)

''' myCNOT
Appends a CNOT to a stim circuit between 
- control qubit (uc, vc, wc)
- target qubit  (ut, vt, wt)
where control and target are tuples'''
def myCNOT(circuit, l, m, control, target, errors: dict):
  uc, vc, wc = control
  ut, vt, wt = target
  kc = convtok(l, m, uc, vc, wc)
  kt = convtok(l, m, ut, vt, wt)
  
  circuit.append("CNOT", [kc, kt])

  p = errors['CNOT'].p
  
  if p > 0:
  
    error_op = errors['CNOT'].op
    circuit.append(error_op, [kc, kt], p)


''' add_A_CNOTs
For A in Hx = [A|B] CNOTs between X check and L data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix A which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CNOTs between each X check qubit (X, v, w) and L data qubit (L, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_A_CNOTs(circuit, jval, code, registers, errors, idle_during, sequential_gates):

  Aij = code.Aij
  l = code.l
  m = code.m
  X = registers.X
  L = registers.L
  qR = registers.qR
  
  for (i, j) in Aij:
      if j == jval:  # i.e. this (i, j) appears in Aij
          for v in range(l): # for each check qubit within a module
              for w in range(m): # for each module
                control = (X, v, w)
                target = (L, (v + i) % l , (w + j) % m) # LEFT qubits as we're doing matrix A in Hx = [A|B]
                myCNOT(circuit, l, m, control, target, errors)
                
              if sequential_gates: # if we are doing one CNOT per timestep (per module) we need to add idling errors to qubits that weren't in the CNOT, namely all the R data qubits and any L data qubits with v' ≠ v ⊕ i

                idle(circuit, qR, idle_during['CNOT']) # idle all the R data qubits
                for w in range(m):
                  for vprime in range(l): # idle L qubits not in CNOT:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, L, vprime, w)
                      idle(circuit, [qubit], idle_during['CNOT'])
                  
                tick(circuit)

          if not sequential_gates: # idle all the R data qubits only after all L data qubits had CNOTs in a single timestep
            idle(circuit, qR, idle_during['CNOT']) # idle the R data qubits
            tick(circuit)

''' add_B_CNOTs
For B in Hx = [A|B], this function appends CNOTs between X check and R data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix B which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CNOTs between each X check qubit (X, v, w) and R data qubit (R, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we run through each v to v ⊕ i (modulo l)'''
def add_B_CNOTs(circuit, jval, code, registers, errors, idle_during, sequential_gates):

  Bij = code.Bij
  l = code.l
  m = code.m
  X = registers.X
  R = registers.R
  qL = registers.qL
  
  for (i, j) in Bij:
      if j == jval:  # i.e. this (i, j) appears in Bij
          for v in range(l):
              for w in range(m):
                control = (X, v, w)
                target = (R, (v + i) % l , (w + j) % m) # RIGHT qubits as we're doing matrix B in Hx = [A|B]
                myCNOT(circuit, l, m, control, target, errors)

              if sequential_gates: # if we are doing one CNOT per timestep (per module) we need to add idling errors to qubits that weren't in the CNOT, namely all the L data qubits and any R data qubits with v' ≠ v ⊕ i

                idle(circuit, qL, idle_during['CNOT']) # idle all the L data qubits
                for w in range(m):
                  for vprime in range(l): # idle R qubits not in CNOT:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, R, vprime, w)
                      idle(circuit, [qubit], idle_during['CNOT'])

                tick(circuit)

          if not sequential_gates: # idle all the L data qubits only after all R data qubits had CNOTs in a single timestep
            idle(circuit, qL, idle_during['CNOT']) # idle the L data qubits
            tick(circuit)




''' myCZ
Appends a CZ to a stim circuit between 
- control qubit (uc, vc, wc)
- target qubit  (ut, vt, wt)
where control and target are tuples'''
def myCZ(circuit, l, m, control, target, errors: dict):
  uc, vc, wc = control
  ut, vt, wt = target
  kc = convtok(l, m, uc, vc, wc)
  kt = convtok(l, m, ut, vt, wt)
  
  circuit.append("CZ", [kc, kt])

  p = errors['CZ'].p

  if p > 0:
    
    error_op = errors['CZ'].op

    circuit.append(error_op, [kc, kt], p)


''' add_BT_CZs
For B^T in Hz = [B^T|A^T], this function appends CZs between Z-check and L data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix B^T which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CZs between each Z check qubit (Z, v, w) and L data qubit (L, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_BT_CZs(circuit, jval, code, registers, errors, idle_during, sequential_gates):

  BTij = code.BTij
  l = code.l
  m = code.m
  Z = registers.Z
  L = registers.L
  qR = registers.qR
  
  for (i, j) in BTij:
      if j == jval:  # i.e. this (i, j) appears in BTij
          for v in range(l): # for each check qubit within a module
              for w in range(m): # for each module
                control = (Z, v, w)
                target = (L, (v + i) % l , (w + j) % m) # LEFT qubits as we're doing matrix BT in Hz = [B^T|A^T]

                myCZ(circuit, l, m, control, target, errors)

                
              if sequential_gates: # if we are doing one CZ per timestep (per module) we need to add idling errors to qubits that weren't in the CZ, namely all the R data qubits and any L data qubits with v' ≠ v ⊕ i

                idle(circuit, qR, idle_during['CZ'])  # idle all the R data qubits
                for w in range(m):
                  for vprime in range(l): # idle L qubits not in CZ:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, L, vprime, w)
                      idle(circuit, [qubit], idle_during['CZ']) 
                tick(circuit)

          if not sequential_gates: # idle all the R data qubits only after all L data qubits had CZs in a single timestep
            idle(circuit, qR, idle_during['CZ'])  # idle the R data qubits
            tick(circuit)




''' add_AT_CZs
For A^T in Hz = [B^T|A^T], this function appends CZs between Z check qubits and R data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix A^T which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CZs between each Z check qubit (Z, v, w) and R data qubit (R, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_AT_CZs(circuit, jval, code, registers, errors, idle_during, sequential_gates):

  ATij = code.ATij
  l = code.l
  m = code.m
  Z = registers.Z
  R = registers.R
  qL = registers.qL
  
  for (i, j) in ATij:
      if j == jval:  # i.e. this (i, j) appears in ATij
          for v in range(l): # for each check qubit within a module
              for w in range(m): # for each module
                control = (Z, v, w)
                target = (R, (v + i) % l , (w + j) % m) # RIGHT qubits as we're doing matrix A^T in Hz = [B^T|A^T]

                myCZ(circuit, l, m, control, target, errors)
                
              if sequential_gates: # if we are doing one CZ per timestep (per module) we need to add idling errors to qubits that weren't in the CZ, namely all the L data qubits and any R data qubits with v' ≠ v ⊕ i 

                idle(circuit, qL, idle_during['CZ'])  # idle all the L data qubits
                for w in range(m):
                  for vprime in range(l): # idle R qubits not in CZ:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, R, vprime, w)
                      idle(circuit, [qubit], idle_during['CZ']) 
                tick(circuit)

          if not sequential_gates: # idle all the L data qubits only after all R data qubits had CZs in a single timestep
            idle(circuit, qL, idle_during['CZ'])  # idle the L data qubits
            tick(circuit)


''' update_shift_probs
The cyclic shift required to align check qubit and data qubit modules can be of variable length depending on the previous and current value of j (recall we are aligning check qubit module M_w with data qubit module M_(w ⊕ j) ). This function takes that length and updates the error rates / probabilities in errors and idle_during dictionaries to correspond to it. It uses an always-stored 'shift_constant' that the error rate is proportional to'''
def update_shift_probs(length_of_shift, errors, idle_during):        
  
  # Get what the shifting error is proportional to:
  p_shift_const = errors['shift_constant']
  
  # Multiply it by the length of the shift (dif. in j values) and update errors dictionary:
  updated_prob = p_shift_const * length_of_shift
  errors['shift'].p = updated_prob

  # Repeat for the idling error:

  p_shift_idle_const = idle_during['shift_constant']

  updated_idle_prob = p_shift_idle_const * length_of_shift

  idle_during['shift'].p = updated_idle_prob



'''apply_cyclic_shifts_and_stab_interactions
Append to a stim circuit for a BB code the reqiured cyclic shifts and two-qubit gates to measure the stabilisers. This is according to Algorithm 2 of Ye ... Delfosse (2 × L array; 2508.01879). Accepts inputs
- circ: the stim circuit to be appended to
- jval_prev: the previous arrangement of modules. If jval_prev = j this implies that check qubit module M^a_0 was aligned with data qubit module M^d_((0 + j) % m)
- check: whether we are performing the X-stabiliser or Z-stabiliser checks
- code: paramaters of the BB code, e.g. Hx = [A|B], Hz = [B^T|A^T] etc.
- registers: indices for the stim circuit of check qubits, data qubits
- errors: what errors and probabilities are on each operation
- sequential (bool): whether the two-qubit gates within a module are applied sequentially (in serial) or in parallel'''
def apply_cyclic_shifts_and_stab_interactions(circ, jval_prev, check, code, registers, errors, idle_during, sequential):

    if check == 'X':
      qC = registers.qX
      theunion = code.Junion
    elif check == 'Z':
      qC = registers.qZ
      theunion = code.JTunion
    else:
        raise ValueError("Parameter 'check' must be either 'X' or 'Z'.")

    l = code.l
    m = code.m 
    qL = registers.qL
    qR = registers.qR



    # Do cyclic shifts to required j-valued modules

    for jval in theunion:

        # # Cyclic shift the check qubits:
        if errors['shift_constant'] > 0:  
  
          length_of_shift = abs((jval % m) - (jval_prev % m))
          update_shift_probs(length_of_shift, errors, idle_during)
          ("updated shift probs")

          apply_shift_error(circ, qC, errors)
          idle(circ, qL + qR, idle_during['shift']) # t_shift) # idle the data qubits 
          tick(circ)

        # Shuttle check qubit modules from racetrack into leg:
        if errors['shuttle'].p > 0:
          apply_shuttle_error(circ, qC, errors)
          idle(circ, qL + qR, idle_during['shuttle'])  # idle the data qubits
          tick(circ)

        # Merge check and data qubit modules Coulomb potentials:
        if errors['merge'].p > 0:
          apply_merge_error(circ, qC + qL + qR, errors)
          tick(circ)

        if check == 'X':
          # Apply CNOTs for X-checks, i.e. Hx = [A|B]
          add_A_CNOTs(circ, jval, code, registers, errors, idle_during, sequential) 
          add_B_CNOTs(circ, jval, code, registers, errors, idle_during, sequential)
        elif check == 'Z':
          # Apply CZs for Z-checks, i.e. Hz = [B^T|A^T]
          add_BT_CZs(circ, jval, code, registers, errors, idle_during, sequential)
          add_AT_CZs(circ, jval, code, registers, errors, idle_during, sequential)

        # Split coulomb potentials of data qubit modules from check qubit modules:
        if errors['split'].p > 0:
          apply_split_error(circ, qC + qL + qR, errors)
          tick(circ)

        # # Shuttle check qubits from leg into racetrack:
        if errors['shuttle'].p > 0:
          apply_shuttle_error(circ, qC, errors)
          idle(circ, qL + qR, idle_during['shuttle']) # idle the data qubits
          tick(circ)

        jval_prev = jval

    return jval_prev



''' make_loop_body
Constructs the stabiliser extraction round of a memory experiment using a BB code. This round or 'loop_body' can be repeated arbitrarily. Returns the loop_body which is a stim circuit. Note this does not initialise data qubits as it is designed to follow an already-constructed round-0 of stabiliser measurements which is slightly different to the repeated rounds (namely round-0 initialises the data qubits and only has detectors on stabilisers in the same basis as the preserved logical state because the other ones are non-deterministic).
- jval_prev: gives the previous arrangement of modules before starting this repeated / looped section. I.e. if jval_prev = j, this implies check module M^a_w was aligned with data module M^d_((w + j) % m)
- code: paramaters of the BB code. Hx = [A|B], Hz = [B^T,A^T] etc.
- errors: error rates and operations for each circuit operation
- memory_basis: 'Z' implies preserve logical 0, 'X' implies preserve logical plus
- exclude_opposite_basis_detectors: if this is True, when preserving logical 0 (memory_basis = 'Z') then there are no detectors placed on the X-stabiliser measurments (though they are still performed). This is useful if the decoder being used is uncorrelated (i.e. treats X and Z detector graphs separately) as it reduces the size of the detector error model to be fed to it by taking out unused nodes.
- reuse_check_qubits: one register of check qubits of size n/2 -- is possible as we're doing X-checks then Z-checks
- sequential: whether or not the two-qubit gates within a leg are sequential or in parallel'''

def make_loop_body(jval_prev, code, errors, idle_during, registers, memory_basis, reuse_check_qubits, sequential, exclude_opposite_basis_detectors = False):

    loop_body = stim.Circuit()

    qX = registers.qX
    qL = registers.qL
    qR = registers.qR
    qZ = registers.qZ  # qZ will equal qX if reuse_check_qubits == True


    # X-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    qC = qX

    # Initialise check qubits
    init('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['RZ']) # t_init) # idle data qubits
    tick(loop_body)

    # Hadamard check qubits to |+⟩
    hadamard(loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['H']) # t_init) # idle data qubits
    tick(loop_body)

    # Do cyclic shifts to required j-valued modules, apply two-qubit gates for stabilisers and return last j position
    jval_prev = apply_cyclic_shifts_and_stab_interactions(loop_body, jval_prev, 'X', code, registers, errors, idle_during, sequential)
    # Alrighty we've done the X-check CNOTs!

    # Hadamard check qubits (which have already been shuttled back to racetrack in apply_cyclic_shifts_and_stab_interactions)
    hadamard(loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['H']) # t_had) # idle data qubits during hadamard
    tick(loop_body)

    # Measure check qubits
    measure('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['MZ']) # t_meas) # idle data qubits
    tick(loop_body)

    # We now place detectors on these check qubit measurements which compare them and the previous round's X-check measurements
    n = code.n
    if (memory_basis == 'Z' and not exclude_opposite_basis_detectors) or memory_basis == 'X':
            # Append X-check stabiliser detectors:
            for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are the X-check measurements, and compares each of them to the measurement performed n measurements before it (this is the same measurement in the preceding round)
                loop_body.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - n)])


    # # Z-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    qC = qZ

    # Initialise Z-check qubits
    init('Z', loop_body, qC, errors) # (note qZ = qX if reuse_check_qubits == True)
    idle(loop_body, qL + qR, idle_during['MZ']) # t_init) # idle data qubits
    tick(loop_body)

    # Hadamard check qubits to |+⟩ and IDLE data qubits:
    hadamard(loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['H']) # t_init) # idle data qubits
    tick(loop_body)

    # Apply required cyclic shifts and CZ interactions for Z-checks:
    jval_prev = apply_cyclic_shifts_and_stab_interactions(loop_body, jval_prev, 'Z', code, registers, errors, idle_during, sequential)

    # Now to hadamard the check qubits (they've already been shuttled back into racetrack)
    hadamard(loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['H'])  # idle data qubits
    tick(loop_body)

    # Now measure check qubits
    measure('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['MZ']) # t_meas) # idle data qubits
    tick(loop_body)

    
    # We now place detectors on these check qubit measurements, comparing them and the previous round's X-check measurements
    n = code.n
    if (memory_basis == 'X' and not exclude_opposite_basis_detectors) or memory_basis == 'Z':
            # Append Z-check stabiliser detectors:
            for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are now Z-check measurements, and compares each of them to the measurement performed n measurements before it (this is the same measurement in the preceding round)
                loop_body.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - n)])

    return loop_body




''' make_circuit
This function makes a stim circuit realising a memory experiment for any bivariate bicycle (BB) code [2308.07915] according to algorithm 2 of Tham et al.'s "qubit modules" paper [2508.01879]. All the X-checks are performed in parallel, then all the Z-checks. The two-qubit gates proceed according to ascending j, the powers of y in the BB code's matrices A = sum(x^i * y^j), B = sum(x^i * y^j), the two matrices which form the code's parity check matrices Hx = [A|B] and Hz = [B^T|A^T]. All of the required information is contained in the input object 'code'.

Inputs are:
    - code
            An object which contains all the BB code's parameters, such as parity check matrices.
    - memory_basis
            Either 'Z' or 'X' to preserve logical 0 or + respectively in all the logical qubits
    - p
            Physical error rate
    - errors
            A dictionary that can be made with errors = default_errors(p) which has an error operation and corresonding probability for each operation in the circuit.
    - idle_during
            A dictionary that can be made with idle_during = default_idle_errors(p) which has an error operation and corresponding probability for qubits that are idling while other qubits are undergoing each operation in the circuit
    - num_syndrome_extraction_cycles
            How many rounds of stabiliser measurements to perform (including the first round which just encodes the logical state)
    - sequential_gates = True
            Whether or not the two-qubit gates between qubits in aligned modules are sequential or in parallel.
    - exclude_opposite_basis_detectors = False
            If this is True, when preserving logical 0 (+) then there are no detectors placed on the X-(Z-) stabiliser measurements (though they are still performed). This is useful if the decoder being used is uncorrelated (i.e. treats X and Z detector graphs separately) as it reduces the size of the detector error model to be fed to it (and thus increases the speed of the simulations) by removing unused detectors.
    - reuse_check_qubits = True
            If true we have one check register of size l * m rather than two, which is reused for X-checks then Z-checks. This is advised and possible because the algorithm this code implements (algorithm 2 of 2508.01879) does X-checks then Z-checks'''
def make_circuit(
  code,
  memory_basis,
  p,
  num_syndrome_extraction_cycles,
  errors = None,
  idle_during = None,
  sequential_gates = True,
  exclude_opposite_basis_detectors = False,
  reuse_check_qubits = True):


    if idle_during == None:
      idle_during = default_idle_errors(p)

    if errors == None:
      errors = default_errors(p)


    circ = stim.Circuit()


    registers = make_registers(code.l, code.m, reuse_check_qubits = reuse_check_qubits)
    qX = registers.qX
    qL = registers.qL
    qR = registers.qR
    qZ = registers.qZ  # qZ will equal qX if reuse_check_qubits == True



    ## Round 0:
    # - wait's a time step to initalise data qubits (we are assuming can't prepare |+⟩ state directly but need RZ then H)
    # - only puts detectors on X (Z) -checks if preparing logical |+⟩ (|0⟩)


    # X-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    qC = qX

    # Initialise X-check qubits
    init('Z', circ, qC, errors)
    tick(circ)

    # Hadamard check qubits to |+⟩ and initalise data qubits to +1 eigenstate of 'memory_basis'
    hadamard(circ, qC, errors)
    init(memory_basis, circ, qL + qR, errors) 
    tick(circ)

    Junion = code.Junion
    jval_0 = Junion[0] # we assume the starting arrangement of the modules is M^a_w with M^d_((w + j) % m), i.e. no cyclic shift errors initially

    # Do cyclic shifts to required j-valued modules (and return last j position)
    jval_prev = apply_cyclic_shifts_and_stab_interactions(circ, jval_0, 'X', code, registers, errors, idle_during, sequential_gates)
    # Alrighty we've done the X-check CNOTs!


    # Now to hadamard the check qubits (they've already been shuttled back into racetrack in apply_cyclic... function)
    hadamard(circ, qC, errors)
    idle(circ, qL + qR, idle_during['H']) # t_had)
    tick(circ)

    # Now measure the check qubits
    measure('Z', circ, qC, errors)
    idle(circ, qL + qR, idle_during['MZ']) # t_meas)
    tick(circ)

    # If preserving logical plus we put detectors on these measurements in the first round:
    n = code.n
    if memory_basis == 'X':
        for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
            circ.append("DETECTOR", [stim.target_rec(-i)])


    # # Z-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    qC = qZ

    # Initialise Z-check qubits
    init('Z', circ, qC, errors) # (note qZ = qX if reuse_check_qubits == True)
    idle(circ, qL + qR, idle_during['MZ']) # t_init, noise) # idle data qubits
    tick(circ)

    # Hadamard check qubits to |+⟩ and IDLE data qubits:
    hadamard(circ, qC, errors)
    idle(circ, qL + qR, idle_during['H']) # t_init, noise) # idle data qubits
    tick(circ)

    # Apply required cyclic shifts and CZ interactions for Z-checks:
    jval_prev = apply_cyclic_shifts_and_stab_interactions(circ, jval_prev, 'Z', code, registers, errors, idle_during, sequential_gates)

    # Now to hadamard the check qubits (they've already been shuttled back into racetrack)
    hadamard(circ, qC, errors)
    idle(circ, qL + qR, idle_during['H'])
    tick(circ)

    # Now measure check qubits
    measure('Z', circ, qC, errors)
    idle(circ, qL + qR, idle_during['MZ']) # t_meas)
    tick(circ)

    # If preserving logical zero we put detectors on these (Z) check measurements in the first round:
    n = code.n
    if memory_basis == 'Z':
        for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
            circ.append("DETECTOR", [stim.target_rec(-i)])

    ## Make repeated / looped stabiliser measurement rounds:
    loop_body = make_loop_body(jval_prev, code, errors, idle_during, registers, memory_basis, reuse_check_qubits, sequential_gates, exclude_opposite_basis_detectors)
    
    # Append loop_body to circuit:
    circ = circ + (num_syndrome_extraction_cycles - 1) * loop_body

    # Final measurement of all data qubits:
    measure(memory_basis, circ, qL + qR, errors)


    ### Add final detectors:
    add_final_detectors(circ, code, memory_basis)


    # Add logical observables (the Lx's or Lz's if mem X or Z):
    add_logical_observables(circ, code.n, code.Lx, code.Lz, memory_basis)


    detecting_regions = circ.detecting_regions() # a test to see that it has valid detecting regions 

    # Save circuit:
    # circ.to_file(f"../circuits/nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},b={memory_basis},noise={noise},r={num_syndrome_extraction_cycles},code=BB,l={l},m={m},A='{''.join(str(x) + str(y) for x, y in Aij)}',B='{''.join(str(x) + str(y) for x, y in Bij)}'.stim")

    # # print(circ.to_crumble_url())
    # svg = str(circ.without_noise().diagram("timeline-svg"))
    # with open("output.svg", "w", encoding="utf-8") as f: f.write(svg)

    return circ