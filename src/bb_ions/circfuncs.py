''' circfuncs
Given the parity check matrices (constructed using bbfuncs.py) and paramaters & logical operators (found using bbparams.py) of a Bicycle Bivariate [2308.07915] code, these functions are for constructing a stim circuit that realises a memory experiment using the BB code.
I.e. prepare logical |0⟩ or |+⟩ in all the logical qubits of the BB code, run multiple rounds of stabiliser measurements, measure all the data qubits.
Measurements results from stabilisers and data qubits will then be given to a decoder to see if it correctly predicts what the logical states were.

Note Stim api reference: https://github.com/quantumlib/Stim/wiki/Stim-v1.9-Python-API-Reference'''

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




''' old_make_registers
  Makes lists of qubit indices dividing 2n qubits evenly into qA, qB, qC, qD.
  This is most commonly used for a BB code, for example qX, qL, qR, qZ, where qX is the X-check syndrome qubits, qL the 'left' data qubits (appear in left-hand side of Hx; acted on by matrix A), qR the 'right' data qubits (appear in right-hand side of Hx; acted on by matrix B) and qZ, the Z-check syndrome qubits.

  qA: qubits 0 to n/2 - 1
  qB: qubits n/2 to n - 1
  qC: qubits n to 3n/2 - 1 
  qD: qubits 3n/2 to 2n - 1'''
def old_make_registers(n):
  qA = list(range( 0 , n//2))
  qB = list(range( n//2 , n))
  qC = list(range( n , 3*n//2))
  qD = list(range( 3*n//2 , 2*n))

  return qA, qB, qC, qD
  


''' initX
Sets qubits in the list 'register' to the plus state |+⟩.
Adds a reset error (set to |-⟩ ) with probability p'''
def initX(circuit, register, p = 0):
  circuit.append("RX", register)
  if p > 0:
    circuit.append("Z_ERROR", register, p)

''' initZ
Sets qubits in the list 'register' to the zero state |0⟩.
Adds a reset error (set to |1⟩ ) with probability p'''
def initZ(circuit, register, p = 0):
  circuit.append("R", register)
  if p > 0:
    circuit.append("X_ERROR", register, p)

''' init
Sets qubits in the list 'register' to
- |0⟩ if basis == 'Z'
- |+⟩ if basis == 'X'
Also adds a reset error (set to |1⟩ or |-⟩ ) with probability p'''
def init(basis, circuit, register, t_init = 0):
  circuit.append(f"R{basis}", register)
  
  p = p_init(t_init)

  if p > 0:
    if basis == 'Z':
      error = 'X_ERROR'
    elif basis == 'X':
      error = 'Z_ERROR'
    
    circuit.append(error, register, p)



''' hadamard
Appends a hadamard gate to the stim circuit on qubit(s) specified in 'register'.
After the gate it adds a depolarising noise of strength p (i.e. an error will occur with prob p. Given it occurs, pick one of X, Y or Z at random)'''
def hadamard(circuit, register, t_had = 0):
  circuit.append("H", register)
  p = p_had(t_had)
  if p > 0:
    circuit.append("DEPOLARIZE1", register, p)


''' measure
Appends a Z-basis measurement onto the qubits specified by 'register' on an input stim circuit 'circuit'. Probability of measurement error given by p_meas(t_meas)'''
def measure(circuit, register, t_meas = 0):
  p = p_meas(t_meas)
  circuit.append("MZ", register, p)

''' tick
Appends a 'TICK' annotation to an input stim circuit, indicating the end of a time-step. '''
def tick(circuit):
  circuit.append("TICK")


''' add_final_detectors
After measuring all the data qubits in the X-basis (memory X) or Z-basis (memory Z) we want to check that each check qubit
has correctly reported the parity of the data qubits it was supposed to have measured.
Consequently, we add detectors that include each check qubit's parity multiplied with the parity of all the data qubits it checked.
Note that in mem. Z we only add these detectors to the Z-checks because the Z-measurements don't commute with X-stabiliser
measurements so we cannot check the parity of X-stabilisers from the data qubits that have been measured in Z
(whatever the parity may have been for an X-stab measurement, e.g. -1 for XXXX if |0000⟩ - |1111⟩ , once we collapse to, say, |1111⟩ there is
no more X-parity check).
Conversely for memory X we only add these detectors to the X-checks.

Optional further reading on inner workings of the function: figuring out each stim.target_rec[]:

Most recent rec (rec[-1]) is the n-th data qubit. It was the last measured. This is the last data qubit in qR. Zeroth data qubit, also the
zeroth data qubit in qL, is rec[-n].

So rec:
-1 to -n is data qubit measurements
-(n + 1) to - (n + n//2) is qZ measurements
-(n + n//2 + 1) to -2n is qX measurements

For each check qubit we want the detector to include its parity multiplied with the parity of all the data qubits it checked. These are contained
in the lists x_check_qubits, z_check_qubits (except starting from zero. map to -1 being the n-th qubit and -n being the zeroth qubit, i.e. minus
n from whatever is in the lists)
'''
def add_final_detectors(circuit, n, ones, memory):

  x_check_qubits, z_check_qubits = get_stabiliser_qubit_indices_BB5(ones)

  if memory == 'X': # add detectors to X-checks
    for i in reversed(list(range(n//2))):
      circuit.append(
          "DETECTOR",
          [ stim.target_rec(-2 * n + i) ] # the X-check, starting from the first X-check, i.e. qX[0], and adding 1 each time to get to the most recent X-check
          +
          [ stim.target_rec(- j - n) for j in x_check_qubits[i] ] # the data qubits it checks the parity of. Minus n off each 0 to n - 1 data qubit indices in x_check_qubits so the zeroth qubit is rec[-n]
          #, [, y, time] # coordinates
      )

  elif memory == 'Z': # add detectors to Z-checks
    for i in reversed(list(range(n//2))):
      circuit.append(
          "DETECTOR",
          [stim.target_rec(- (n + n//2) + i)] # the Z-check, starting from the first Z-check (NOT the most recent)
          +
          [ stim.target_rec(- j - n) for j in z_check_qubits[i] ] # all the data qubits it checks. Minus n off each 0 to n - 1 data qubit index
      )


''' add_logical_observables
The circuit ends with a parity check of the logical operators / observables. In a surface code, for example, XL and ZL are just vertical or horizontal chains of X's and Z's across the lattice.
In a BB code they are also chains of X's and Z's but on specific qubits, contained in the arrays Lx and Lz. For each logical qubit Lx and Lz contain a pair of anti-commuting logical operators.
These commute with the logical operators of other logical qubits. This function adds the Lx operators as observable if we are in memory X (preserving an eigenstate of Lx's, i.e. |+⟩_L)
or the Lz operators if we are in memory Z.

Optional note on inner workings of function:
The data qubit measurements were the most recent measurements in the circuit so are from rec[-1] (the n-th data qubit) to rec[-n] (the zeroth data qubit). So if Lx has a 1 on qubit 3,
for example, this means it needs rec[3 - n]. Just minus n from the index in Lx.

'''
def add_logical_observables(circuit, n, Lx, Lz, memory):

  num_logical_ops = Lx.shape[0]

  L = Lx if memory == 'X' else Lz

  indices = get_nonzero_indices(L) # instead of L being 1's and 0's (like a parity check matrix) just make it a list of the indices of the 1's

  for i in range(num_logical_ops): # for each logical op

    recordings = (indices[i] - n).astype(int) # the measurements -- 'inner workings' note above

    circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(r) for r in recordings], 0.0)

''' myCNOT
Appends a CNOT to a stim circuit between 
- control qubit (uc, vc, wc)
- target qubit  (ut, vt, wt)
where control and target are tuples'''
def myCNOT(circuit, l, m, control, target, t_cnot):
  uc, vc, wc = control
  ut, vt, wt = target
  kc = convtok(l, m, uc, vc, wc)
  kt = convtok(l, m, ut, vt, wt)

  p = p_cnot(t_cnot)
  
  circuit.append("CNOT", [kc, kt])

  if p > 0:
    circuit.append("DEPOLARIZE2", [kc, kt], p)


''' add_A_CNOTs
For A in Hx = [A|B] CNOTs between X check and L data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix A which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CNOTs between each X check qubit (X, v, w) and L data qubit (L, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_A_CNOTs(circuit, jval, code, registers, t_cnot, sequential_gates):

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
                myCNOT(circuit, l, m, control, target, t_cnot)
                
              if sequential_gates: # if we are doing one CNOT per timestep (per module) we need to add idling errors to qubits that weren't in the CNOT, namely all the R data qubits and any L data qubits with v' ≠ v ⊕ i

                idle(circuit, qR, t_cnot) # idle all the R data qubits
                for w in range(m):
                  for vprime in range(l): # idle L qubits not in CNOT:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, L, vprime, w)
                      idle(circuit, [qubit], t_cnot)
                  
                tick(circuit)

          if not sequential_gates: # idle all the R data qubits only after all L data qubits had CNOTs in a single timestep
            idle(circuit, qR, t_cnot) # idle the R data qubits
            tick(circuit)

''' add_B_CNOTs
For B in Hx = [A|B], this function appends CNOTs between X check and R data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix B which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CNOTs between each X check qubit (X, v, w) and R data qubit (R, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we run through each v to v ⊕ i (modulo l)'''
def add_B_CNOTs(circuit, jval, code, registers, t_cnot, sequential_gates):

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
                myCNOT(circuit, l, m, control, target, t_cnot)

              if sequential_gates: # if we are doing one CNOT per timestep (per module) we need to add idling errors to qubits that weren't in the CNOT, namely all the L data qubits and any R data qubits with v' ≠ v ⊕ i

                idle(circuit, qL, t_cnot) # idle all the L data qubits
                for w in range(m):
                  for vprime in range(l): # idle R qubits not in CNOT:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, R, vprime, w)
                      idle(circuit, [qubit], t_cnot)

                tick(circuit)

          if not sequential_gates: # idle all the L data qubits only after all R data qubits had CNOTs in a single timestep
            idle(circuit, qL, t_cnot) # idle the L data qubits
            tick(circuit)




''' myCZ
Appends a CZ to a stim circuit between 
- control qubit (uc, vc, wc)
- target qubit  (ut, vt, wt)
where control and target are tuples'''
def myCZ(circuit, l, m, control, target, t_cz):
  uc, vc, wc = control
  ut, vt, wt = target
  kc = convtok(l, m, uc, vc, wc)
  kt = convtok(l, m, ut, vt, wt)

  p = p_cz(t_cz)
  
  circuit.append("CZ", [kc, kt])

  if p > 0:
    circuit.append("DEPOLARIZE2", [kc, kt], p)


''' add_BT_CZs
For B^T in Hz = [B^T|A^T], this function appends CZs between Z-check and L data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix B^T which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CZs between each Z check qubit (Z, v, w) and L data qubit (L, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_BT_CZs(circuit, jval, code, registers, t_cz, sequential_gates):

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

                myCZ(circuit, l, m, control, target, t_cz)

                
              if sequential_gates: # if we are doing one CZ per timestep (per module) we need to add idling errors to qubits that weren't in the CZ, namely all the R data qubits and any L data qubits with v' ≠ v ⊕ i

                idle(circuit, qR, t_cz) # idle all the R data qubits
                for w in range(m):
                  for vprime in range(l): # idle L qubits not in CZ:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, L, vprime, w)
                      idle(circuit, [qubit], t_cz)
                tick(circuit)

          if not sequential_gates: # idle all the R data qubits only after all L data qubits had CZs in a single timestep
            idle(circuit, qR, t_cz) # idle the R data qubits
            tick(circuit)




''' add_AT_CZs
For A^T in Hz = [B^T|A^T], this function appends CZs between Z check qubits and R data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix A^T which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies CZs between each Z check qubit (Z, v, w) and R data qubit (R, v ⊕ i, w ⊕ j) for one value of j. The required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), now within modules we do each v to v ⊕ i (modulo l)'''
def add_AT_CZs(circuit, jval, code, registers, t_cz, sequential_gates):

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

                myCZ(circuit, l, m, control, target, t_cz)
                
              if sequential_gates: # if we are doing one CZ per timestep (per module) we need to add idling errors to qubits that weren't in the CZ, namely all the L data qubits and any R data qubits with v' ≠ v ⊕ i 

                idle(circuit, qL, t_cz) # idle all the L data qubits
                for w in range(m):
                  for vprime in range(l): # idle R qubits not in CZ:
                    if vprime != (v + i) % l:
                      qubit = convtok(l, m, R, vprime, w)
                      idle(circuit, [qubit], t_cz)
                tick(circuit)

          if not sequential_gates: # idle all the L data qubits only after all R data qubits had CZs in a single timestep
            idle(circuit, qL, t_cz) # idle the L data qubits
            tick(circuit)