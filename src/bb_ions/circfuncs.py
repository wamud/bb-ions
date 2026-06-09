''' circfuncs
Given the parity check matrices (constructed using bbfuncs.py) and paramaters & logical operators (found using bbparams.py) of a Bicycle Bivariate [2308.07915] code, these functions (culminating in make_BB_circuit at the bottom of this file) are for constructing a stim circuit that realises, with appropriate detectors annotated, a memory experiment using the BB code, i.e. it prepares logical |0⟩ or |+⟩ in all the logical qubits of the BB code, runs multiple rounds of stabiliser measurements, then measures all the data qubits.
Note Stim api reference: https://github.com/quantumlib/Stim/wiki/Stim-v1.9-Python-API-Reference'''

import stim
from .kfuncs import *
from .noisefuncs import *
import math
from typing import Any

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


''' order_of_ops
The syndrome extraction circuit we are following is as per Algorithm 2 in Edwin Tham et al.'s paper [2508.01879]. The polynomials in the BB code are A and B and made of sums of terms x^i⋅y^j. We go through each value of j and then apply any values of i that appear. This function prints that out.'''
def order_of_ops(code):

    l = code.l
    m = code.m 

    print('X-checks')
    for j in code.Junion:
        print(f"j = {j}")
        print("  A")
        for (ival, jval) in code.Aij:
            
            if jval == j:
                print(f"  i = {ival}")
        
        print("  B")
        for (ival, jval) in code.Bij:
            
            if jval == j:
                print(f"  i = {ival}")
    
    print('\nZ-checks')
    for j in code.JTunion:
        print(f"j = {j%m}")
        print("  AT")
        for (ival, jval) in code.ATij:
            
            if jval == j:
                print(f"  i = {ival%l}")
        
        print("  BT")
        for (ival, jval) in code.BTij:
            
            if jval == j:
                print(f"  i = {ival%l}")


''' make_registers
Makes lists of qubit indices dividing n = 2lm qubits evenly into qA, qB, qC, qD.
  qA: qubits 0 to n/2 - 1
  qB: qubits n/2 to n - 1
  qC: qubits n to 3n/2 - 1 
  qD: qubits 3n/2 to 2n - 1'''
def make_registers(l, m, reuse_check_qubits = False):
  
  
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


''' add_qubit_coordinates
Adds coordinates to the qubits for use in circuit diagrams and importing to crumble.
Note that (0,0) is the top left of the diagram, and increasing row or column moves down or right respectively.
NOTE: in Stim it is QUBIT_COORD(COLUMN, ROW) 
We will place qubits in blocks of (reading clockwise from top left) X-check, L-data, Z-check, R-data.
Implying that the zeroth X-check qubit at (0,0), zeroth L-data at (0,1), zeroth R-data at (1,0) and zeroth Z-check at (1,1).
We will do l rows by m columns of each of these blocks.
If reuse check qubits is true we will simply not fill in Z-check qubits
Workings:
To convert indices k into rows of length m (i.e. m columns per row) it is:
row = floor(k/m) , column = k mod m
We then additionally need to place each qubit in the appropritate spot in its block so will have some factors of 2'''
def add_qubit_coordinates(circ, code, registers, reuse_check_qubits):
  
  qX = registers.qX
  qL = registers.qL
  qR = registers.qR
  qZ = registers.qZ

  l = code.l
  m = code.m
  n = code.n

  
  # qX: top left of each block:
  for k in range(n//2):
    qubit = qX[k]
    v, w = convtorowcol(m, k)
    x_coord = 2 * w
    y_coord = 2 * v
    # x_coord =  2 * math.floor(k/m)
    # y_coord =  2 * (k % m)
    circ.append("QUBIT_COORDS", qubit, [x_coord, y_coord])
  
  # qL: top right of each block.
  # Will be same row as qX, then to the right by 1 for each column.
  for k in range(n//2):
    qubit = qL[k]
    v, w = convtorowcol(m, k)
    x_coord = 2 * w + 1
    y_coord = 2 * v
    # x_coord =    2 * math.floor(k/m)  
    # y_coord =  (2 * (k % m)) + 1
    circ.append("QUBIT_COORDS", qubit, [x_coord,y_coord])

  # qR: bottom left of each block.
  # Same x_coord of qX, one higher row
  for k in range(n//2):
    qubit = qR[k]
    v, w = convtorowcol(m, k)
    x_coord = 2 * w
    y_coord = 2 * v + 1
    # x_coord =    (2 * math.floor(k/m)) + 1
    # y_coord =  2 * (k % m)
    circ.append("QUBIT_COORDS", qubit, [x_coord, y_coord])

  # qZ: if reusing one set of check qubits we don't have separate qZ and will just draw qX. If NOT reusing then qZ is bottom right of each block.
  # Same as L-data, but add 1 to each of its rows.
  if reuse_check_qubits == False:
    for k in range(n//2):
      qubit = qZ[k]
      v, w = convtorowcol(m, k)
      x_coord = 2 * w + 1
      y_coord = 2 * v + 1
      # x_coord =    (2 * math.floor(k/m)) + 1
      # y_coord =  (2 * (k % m)) + 1
      circ.append("QUBIT_COORDS", qubit, [x_coord, y_coord])

  return 

''' init
Sets qubits in the list 'register' to
- |0⟩ if basis == 'Z'
- |+⟩ if basis == 'X'
Also adds a depolarizing error (in alignment with lonchain paper 2503.2207) with probability p if 'longchain' in the noise parameter, otherwise the more usual reset error which sets to orthog. eigenstate'''
def init(basis, circuit, register, errors: dict):
  
  reset = f"R{basis}"
  circuit.append(reset, register) # append RZ or RX
  
  if LEAKAGE:    # appending leakage AFTER the reset (wouldn't do anything if appended before)
    add_relax_then_leak(reset, circuit, register, errors) # If errors = helios_errors(p) then this won't actually add any leakage / relax because p_leak = p_relax = 0 for reset gates in ions because "Typically, ions are initialized using optical pumping techniques which do not result in leakage (https://doi.org/10.1103/PhysRevA.100.032325)"



  p = errors[reset].p

  if p > 0:
    error_op = errors[reset].op
    circuit.append(error_op, register, p)



''' hadamard
Appends a hadamard gate to the stim circuit on qubit(s) specified in 'register'.
After the gate it adds a depolarising noise of strength p (i.e. an error will occur with prob p. Given it occurs, pick one of X, Y or Z at random)'''
def hadamard(circuit, register, errors: dict):

  if LEAKAGE:
    add_relax_then_leak("H", circuit, register, errors)

  circuit.append("H", register)

  
  p = errors['H'].p
  
  if p > 0:
    circuit.append(errors['H'].op, register, p)


''' measure
Appends a 'basis'-basis measurement onto the qubits specified by 'register' on an input stim circuit 'circuit'. Probability of measurement error given by p_meas(t_meas)'''
def measure(basis: str, circuit, register, errors: dict):
  
  measure_string = f"M{basis}"
  
  p = errors[measure_string].p

  # Append leakage:
  if LEAKAGE:
    add_relax_then_leak(measure_string, circuit, register, errors)
  
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
We optionally also appended a leakage herald to each of these measurements, which would double the number of check-qubit measurements from n to 2n.

if leakage_heralds == False:
  
  The most recent rec (rec[-1]) is the n-th data qubit. It was the last measured. This is the last data qubit in qR. The zeroth data qubit, also the zeroth data qubit in qL, is rec[-n]. Before this is qX then qZ check qubits.
  
  So rec[ ]:
    -1 to -n is data qubit measurements
    -(n + 1) to -(n + n//2) is qZ measurements
    -(n + n//2 + 1) to -2n is qX measurements

elif leakage_heralds == True:

    rec[ ]:
      -1 to -n                        # n leakage heralds on the n data qubit measurements
      -(n + 1) to -2n                 # n data qubit measurements
      -(2n + 1) to -(2n + n//2)       # n/2 leakage heralds on the n/2 qZ measurements
      -(2n + n//2 + 1) to -3n         # n/2 qZ measurements
      -(3n + 1) to (-3n + n//2)       # n/2 leakage heralds on the n/2 qX measurements
      -(3n + n//2 + 1) to -4n         # n/2 qX measurements


For each check qubit we want the detector to include its parity multiplied with the parity of all the data qubits it checked. These are contained in Hz and Hz, but note that the j-th data qubit will be rec[j - n] without leakage heralds or rec[j - 2n] with leakage heralds (so the 0-th data qubit was measured n or 2n measurements ago, the last data qubit, i.e. the (n - 1)th data qubit, was measured one / (1 + n) measurement(s) ago etc.'''
def add_final_detectors(circ, code, memory_basis):

  Hx = code.Hx
  Hz = code.Hz
  n = code.n

  offset_to_first_data_qubit_measurement_recording = 2 * n if LEAKAGE_HERALDS else n

  if memory_basis == 'Z':
      
      for k in list(range(n//2)): # loop over Z-check qubits
          
          this_check_qubits_data_qubits = np.nonzero(Hz[k])[0] # extract the 1 positions (data qubit indices) in this check qubit's row of Hz

          this_z_check_qubit_rec = -(n + n//2) + k if not LEAKAGE_HERALDS else -3 * n + k     # See above for explanation. Also note we start from the most negative recording and approach the most recent (-1) by adding k each iteration

          circ.append("DETECTOR", 
            [stim.target_rec(this_z_check_qubit_rec)] # This Z-check qubit
            +
            [stim.target_rec(j - offset_to_first_data_qubit_measurement_recording) for j in this_check_qubits_data_qubits] # Its data qubits
            # Note it is j - n (without leakage heralds) because if j is the last data qubit, j = n - 1, we need stim.target_rec(-1) (most recent measurement)
            , [k, k]
            ) # putting at coordinates k,k to make slighty more compact timeline-svg. Not sure what it's doing for other diagrams.

  elif memory_basis == 'X':
      
      for k in list(range(n//2)):

          this_check_qubits_data_qubits = np.nonzero(Hx[k])[0]

          this_x_check_qubit_rec = -2 * n + k if not LEAKAGE_HERALDS else -4 * n + k

          circ.append("DETECTOR", 
            [stim.target_rec(this_x_check_qubit_rec)] 
            +
            [stim.target_rec(j - offset_to_first_data_qubit_measurement_recording) for j in this_check_qubits_data_qubits] # Its data qubits
            ,[k, k]
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

  indices = get_nonzero_indices(L) 

  offset = n if LEAKAGE_HERALDS else 0 # (account for the n data qubit measurements having leakage heralds on them or not)

  for i in range(num_logical_ops): # for each logical qubit

    recordings = (indices[i] - n - offset).astype(int) # the measurements -- 'inner workings' note above

    circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(r) for r in recordings], i)




''' myCP
My 'controlled-Pauli' function. Adds a CX or CZ to a stim circuit between 
- control qubit (uc, vc, wc)
- target qubit  (ut, vt, wt)
where control and target are tuples (see make_registers for index to tuple)'''
def myCP(circuit, gate, l, m, control, target, errors: dict):

  
  
  uc, vc, wc = control
  ut, vt, wt = target
  kc = convtok(l, m, uc, vc, wc)
  kt = convtok(l, m, ut, vt, wt)

  if LEAKAGE_REPUMPING: # Going to repump BEFORE each two-qubit gate (will account for leakage from transport)
    if REPUMPING_CYCLES != 0:
      # cycles = 4 # num_repumping_cycles 
      p_relax = 1 - 1/(3 ** REPUMPING_CYCLES) # based on Honeywell (now Quantinuum after merger) paper 10.1103/PhysRevLett.124.170501
      p_error = REPUMPING_CYCLES * 2e-5       #    ditto
      circuit.append("RELAX", [kc, kt], p_relax)
      circuit.append("DEPOLARIZE1", [kc, kt], p_error)
  
  if LEAKAGE:
    add_relax_then_leak(gate, circuit, [kc, kt], errors)


  circuit.append(gate, [kc, kt])

  p = errors[gate].p

  if p != 0: # using != to account for when p is a list of probs for PAULI_CHANNEL_2
    
    error_op = errors[gate].op

    circuit.append(error_op, [kc, kt], p)

  


  


'''add_2q_gates
Works when the qubit registers are updating due to SWAP gates. Apart from that it's the same as add_2q_gates.

Notes from add_2q_gates function: 
For matrix M (A, B, AT or BT) in Hx = [A|B] or Hz = [B^T|A^T], this adds the chosen two-qubit gate between check qubits and data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix M which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies 2q gates between each check qubit (v, w) and data qubit (v ⊕ i, w ⊕ j) for one value of j. A (B) means X-checks to L(R)-data, BT (AT) means Z-check qubits to L(R)-data qubits. 
We are imagingin that alignging the required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), so now within modules we just do each v to v ⊕ i (modulo l).
If there are sequential gates then a shuttling time is added between each 2q gate.
Addition: swap_ctrl_target variable. If swap == True this swaps the control and target of the two-qubit gate. Useful for doing a leakage reduction circuit using swap gates (SWAP-LRC) as this requires the addition of only one extra CNOT in the opposite direction to the final one (when working only in CNOT gates) or requires the addition of a CZ gate sandwiched by Hadamards (Hadamards added outside this function)'''
def add_2q_gates(gate, matrix, circuit, jval, code, registers, errors, idle_during, sequential_gates, swap_ctrl_target = False):

  if gate == 'CX':
    gate = 'CNOT'
  
  l = code.l
  m = code.m
  # Usually 0, 1, 2, 3:
  X = registers.X
  L = registers.L
  R = registers.R
  Z = registers.Z
  # Qubits contained therein:
  qX = registers.qX
  qL = registers.qL
  qR = registers.qR
  qZ = registers.qZ

  # if swap_ctrl_target == True and gate == 'CZ':
  #   raise ValueError("Swapping control and target only significant with CX gate")


  Mij = getattr(code, f"{matrix}ij")


  # Hx = [A|B], Hz = [BT|AT]
  
  if matrix == 'A' or matrix == 'BT':
    D = L
    qD = qL
    qD_idle = qR

  if matrix == 'B' or matrix == 'AT':
    D = R
    qD = qR
    qD_idle = qL

  if matrix == 'AT' or matrix == 'BT':
    qC = qZ
  if matrix == 'A' or matrix == 'B':
    qC = qX 
  
  

  for (i, j) in Mij: 
      if j == jval: 

        idx_list = []
        a_list = []
        
        for idx in range(l * m):
          
          if matrix == 'A' or matrix == 'B': # Hx = [A|B]
            
            ## Updated code which gets updated qubit indices:
            check_qubit = convtouvw(l, m, qC[idx])

            # Let's say we're dealing with index 0. If (i, j) = (0, 0) we will just interact with qD[0]. If (i,j) = (2, 1) we need to interact with qD[(2, 1)]
            v, w = convtorowcol(m, idx)
            # Shift the index to interact with the (v+i, w+j)-th data qubit
            a = conv_vw_to_k((v + i) % l, (w + j) % m, m)
            a_list.append(a)
            data_qubit = convtouvw(l, m, qD[a])

            control = check_qubit
            target = data_qubit


          if matrix == 'BT' or matrix == 'AT': # Hz = [BT|AT] 

            check_qubit = convtouvw(l, m, qC[idx]) # will be target

            # Let's say we're dealing with index 0. If (i, j) = (0, 0) we will just interact with qD[0]. If (i,j) = (2, 1) we need to interact with qD[(2, 1)]
            v, w = convtorowcol(m, idx)
            # Shift the index to interact with the (v+i, w+j)-th data qubit
            a = conv_vw_to_k((v + i) % l, (w + j) % m, m)
            a_list.append(a)
            data_qubit = convtouvw(l, m, qD[a])

            control = data_qubit
            target = check_qubit

          if swap_ctrl_target == True:
            control_temp = control
            control = target
            target = control_temp

          # 
          # print(f"Control = {control}")
          # print(f"Target = {target}")

          myCP(circuit, gate, l, m, control, target, errors)

          idx_list.append(idx)

          if sequential_gates:
            if (idx + 1) % m == 0: # We can do m gates per time step. Could replace to be more or less but needs to be a multiple of m for the layout. For our layout of m modules we are saying we have applied one 2q gate per module. Assuming sequential gates, i.e. we have m modules and only m operation zones which do one 2q gate per timestep, we must therefore apply idling here to any qubits not in the 2q gate.
              
              idle(circuit, qD_idle, idle_during[gate]) # all the data qubits not in this term at all (e.g. L-data qubits if we're connecting to R-data qubits)

              check_idles = []
              data_idles = []
              
              for g in range(l * m):
                if g not in idx_list:
                  check_idles.append(g)
                if g not in a_list:
                  data_idles.append(g)
              
              check_qubits_to_idle = [qC[entry] for entry in check_idles]
              data_qubits_to_idle = [qD[entry] for entry in data_idles]
              qubits_to_idle = check_qubits_to_idle + data_qubits_to_idle
              
              idle(circuit,qubits_to_idle, idle_during[gate])

                # if g not in idx_list: # for the check qubits, the qubits operated on are simply qC[idx_list], so idle the ones not in idx_list: 
                #   idle(circuit, [qC[g]], idle_during[gate])
                # if g not in a_list: # data qubits
                #   idle(circuit, [qD[g]], idle_during[gate])

              idx_list = []  # reset it 
              a_list = []
              check_idles = []
              data_idles = []

              tick(circuit)


        if not sequential_gates: # idle all data qubits not in this matrix.
          idle(circuit, qD_idle, idle_during[gate])
          tick(circuit)


'''add_2q_gates
The same as add_2q_gates but have an extra if statement to only apply 2q gates for a specific ij value.

For matrix M (A, B, AT or BT) in Hx = [A|B] or Hz = [B^T|A^T], this adds the chosen two-qubit gate between check qubits and data qubits according the the value of j (indicating the modules that are aligned) and the terms in matrix M which have y^j.
As per Algo. 2 & lemma 1 of [2508.01879], this applies 2q gates between each check qubit (v, w) and data qubit (v ⊕ i, w ⊕ j) for one value of j. A (B) means X-checks to L(R)-data, BT (AT) means Z-check qubits to L(R)-data qubits. 
We are imagingin that alignging the required w to w ⊕ j (modulo m) has already been taken care of by aligning modules (simulated by applying required noise), so now within modules we just do each v to v ⊕ i (modulo l).
If there are sequential gates then a shuttling time is added between each 2q gate.
Addition: swap_ctrl_target variable. If swap == True this swaps the control and target of the CNOTs. Useful for doing a leakage reduction circuit using swap gates (SWAP-LRC) as this requires the addition of only one extra CNOT in the opposite direction to the final one (when working only in CNOT gates)'''
def add_2q_gates_for_this_ij_value(gate, matrix, circuit, ival, jval, code, registers, errors, idle_during, sequential_gates, swap_ctrl_target = False):
  
  if gate == 'CX':
    gate = 'CNOT'
  
  l = code.l
  m = code.m
  # Usually 0, 1, 2, 3:
  X = registers.X
  L = registers.L
  R = registers.R
  Z = registers.Z
  # Qubits contained therein:
  qX = registers.qX
  qL = registers.qL
  qR = registers.qR
  qZ = registers.qZ

  # if swap_ctrl_target == True and gate == 'CZ':
  #   raise ValueError("Swapping control and target only significant with CX gate")


  Mij = getattr(code, f"{matrix}ij")


  # Hx = [A|B], Hz = [BT|AT]
  
  if matrix == 'A' or matrix == 'BT':
    D = L
    qD = qL
    qD_idle = qR

  if matrix == 'B' or matrix == 'AT':
    D = R
    qD = qR
    qD_idle = qL

  if matrix == 'A' or matrix == 'B':
    qC = qX
  if matrix == 'AT' or matrix == 'BT':
    qC = qZ

  
  for (i, j) in Mij: # runs over all i,j, applies any i that has this value of j:
      if j == jval: # i.e. this (i, j) appears in Mij, we apply this value of i:
        if i == ival:
          
          idx_list = []
          a_list = []
          
          for idx in range(l * m):
            

            # Target and control set as if doing all CNOT gates. For CZ gates target / control doesn't matter.
            
            if matrix == 'A' or matrix == 'B': # Hx = [A|B]

              ## Updated code which gets updated qubit indices:
              check_qubit = convtouvw(l, m, qC[idx])

              # Let's say we're dealing with index 0. If (i, j) = (0, 0) we will just interact with qD[0]. If (i,j) = (2, 1) we need to interact with qD[(2, 1)]
              v, w = convtorowcol(m, idx)
              # Shift the index to interact with the (v+i, w+j)-th data qubit
              a = conv_vw_to_k((v + i) % l, (w + j) % m, m)
              a_list.append(a)
              data_qubit = convtouvw(l, m, qD[a])

              control = check_qubit
              target = data_qubit


            if matrix == 'BT' or matrix == 'AT': # Hz = [BT|AT] 

              check_qubit = convtouvw(l, m, qC[idx]) # will be target

              # Let's say we're dealing with index 0. If (i, j) = (0, 0) we will just interact with qD[0]. If (i,j) = (2, 1) we need to interact with qD[(2, 1)]
              v, w = convtorowcol(m, idx)
              # Shift the index to interact with the (v+i, w+j)-th data qubit
              a = conv_vw_to_k((v + i) % l, (w + j) % m, m)
              a_list.append(a)
              data_qubit = convtouvw(l, m, qD[a])

              control = data_qubit
              target = check_qubit


            if swap_ctrl_target == True:
              control_temp = control
              control = target
              target = control_temp

            myCP(circuit, gate, l, m, control, target, errors)

            idx_list.append(idx)

            if sequential_gates:
              if (idx + 1) % m == 0: # This is assuming we can do m gates per time step (replace m for more or less gates). For our layout of m modules we are saying we have applied one 2q gate per module. Assuming sequential gates, i.e. we have m modules and only m operation zones which do one 2q gate per timestep, we must therefore apply idling here to any qubits not in the 2q gate.
                idle(circuit, qD_idle, idle_during[gate]) # all the data qubits not in this term at all
                
                for g in range(l * m):
                  if g not in idx_list: # for the check qubits, the qubits operated on are simply qC[idx_list], so idle the ones not in idx_list: (for data qubits will be more subtle)
                    idle(circuit, [qC[g]], idle_during[gate])
                  if g not in a_list:
                    idle(circuit, [qD[g]], idle_during[gate])

                idx_list = []  # reset it 
                a_list = []
                
                tick(circuit)


  if not sequential_gates: # The above for loop only excecutes one particular (i,j) value. This idling on the non-touched data qubits consequently happens outside the for loop (rather than within it as in add_2q_gates function which applies multiple i values - the for loop executes something for every i value)
    idle(circuit, qD_idle, idle_during[gate])
    tick(circuit)



''' update_shift_probs
The cyclic shift required to align check qubit and data qubit modules can be of variable length depending on the previous and current value of j (recall we are aligning check qubit module M_w with data qubit module M_(w ⊕ j) ). This function takes that length and updates the error rates / probabilities in errors and idle_during dictionaries to correspond to it. It uses an always-stored 'shift_prop_to' that the error rate is proportional to.

If the 'shift_prop_to' in the errors and/or idle_during dictionaries is set to None then the 'shift' probability in the errors and/or idle_during dictionaries will not be updated, it remains whatever it was initially set to (for example that Tham, Ye .. Delfosse modules paper has it set to a constant 30p/100)'''
def update_shift_probs(j_dif, errors, idle_during):        
  
  # Get what the shifting error is proportional to -- at the moment this is T2
  shift_prop_to = errors['shift_prop_to'] # T2 time

  if shift_prop_to != None: # E.g. in Tham modules noise the shift constant is none, meaning there is no update depending on the length of the shift -- Tham ... Delfosse assumed constant shift error of 30p/100 regardless of length of shift. 
    
    # # Usual function: multiply it by the length of the shift (dif. in j values) and update errors dictionary:
    # updated_prob =  shift_prop_to * j_dif

    # # New modification: for our architecture, where we just say the shift error is equivalent to idling for that length of time
    # # TO DO: work into making an option able to be specified when make_BB_circuit() is called. For now just manually adding it to make the circuits and start running them.
    
    length_of_shift = LEG_SPACING * j_dif  # LEG_SPACING defined in noisefuncs
    shuttle_speed = 1 # [m/s]
    t = length_of_shift / shuttle_speed  # [s]
    updated_prob = (1 / 2) * (1 - np.exp(- t / shift_prop_to)) # A dephasing noise channel with T2 = shift_prop_to

    errors['shift'].p = updated_prob

  # Repeat for the idling error:

  p_shift_idle_const = idle_during['shift_prop_to']

  if p_shift_idle_const != None:
    
    # # Previous:
    # updated_idle_prob = p_shift_idle_const * j_dif
    
    
    # New modification:
    length_of_shift = LEG_SPACING * j_dif 
    shuttle_speed = 1 # m/s
    t = length_of_shift / shuttle_speed # s
    updated_idle_prob = (1 / 2) * (1 - np.exp(- t / shift_prop_to)) # A dephasing noise channel with T2 = shift_prop_to

    idle_during['shift'].p = updated_idle_prob

'''update_qubit_indices
If implementing a swap-LRC, each check qubit is swapped with the last data qubit it interacts with before the new check qubits are measured. This implies the data qubits are now on the check register, however their indexes are shifted by whatever the last interaction was becuase check qubit (v,w) interacted with data qubit (v+i, w+j) for term x^i⋅y^j.
So in subsequent checks interactions need to between what were the data qubits (now check qubits) with the correct data qubits. 

The best way to do this is simply updating the integers X, L, R, Z and the qubit indices contained in qX, qL, qR, qZ. 
So X, L, R, Z will take on different integer values, then the indices within qX, qL, qR, qZ will change according to the swap gate they went through.
For example, if the last term in the X-checks is in A and is x^2 y^3  then this means check qubit (X, v, w) interacted with and was swapped with (L, (v + i) mod l, (w + j) mod m).

This function does the required update and is used in the function add_swap_lrc

matrix -- the matrix (A,B,AT or BT) that the term was in
i, j -- the i and j of the term x^i⋅y^j 
check -- are we doing'X' or 'Z' stabiliser check'''
def update_qubit_indices(code, registers, matrix, check, i, j, reuse_check_qubits):

  # Get current data qubits:
  if matrix == 'A' or matrix == 'BT':
    D = registers.L
    qD = registers.qL
  elif matrix == 'B' or matrix == 'AT':
    D = registers.R
    qD = registers.qR

  # Get current check qubits:
  if check == 'X':
    C = registers.X
    qC = registers.qX
  elif check == 'Z':
    C = registers.Z
    qC = registers.qZ

  # Get code paramaters
  l = code.l
  m = code.m 
  
  # Save current registers
  temp_qC = qC.copy()
  temp_qD = qD.copy()

  # Start updating qubit registers
  for k in range(l * m):

    # we want to swap qC[k], where k = (v, w), with data qubit qD[k'], where k' = (v + i , w + j)
    v, w = convtorowcol(m, k)
    vprime, wprime = (v + i) % l, (w + j) % m
    kprime = conv_vw_to_k(vprime, wprime, m)
    # update the registers:
    qC[k] = temp_qD[kprime] # replacing check qubit k with data qubit k'
    qD[kprime] = temp_qC[k]      # replacing data qubit k' with check qubit k

  # Update register indices:
  temp_C = C
  C = D
  D = temp_C

  # So we now have an updated qC, qD, C and D. So now we want to replace what is in the registers object for qC, qD, C and D (with the exact C ∈ {X, Z} and D ∈ {L, R} depending on what matrix and check we did) with these updated values:
  if matrix == 'A' or matrix == 'BT':
    registers.L = D
    registers.qL = qD
  elif matrix == 'B' or matrix == 'AT':
    registers.R = D
    registers.qR = qD
  
  if check == 'X':
    registers.X = C
    registers.qX = qC
    if reuse_check_qubits == True:
      registers.Z = C
      registers.qZ = qC
  
  elif check == 'Z':
    registers.Z = C
    registers.qZ = qC
    if reuse_check_qubits == True:
      registers.X = C
      registers.qX = qC




def add_hadamards_before_swap_CZ(matrix, thegate, check, registers, circ, errors, idle_during):
  if thegate != 'CZ':
    raise ValueError("This is written to go around a CZ gate")
  # to do a reversed CNOT we sandwich both ends of the CZ by Hadamards. This is on the target in order to make the gate a CNOT. It is also on the control in order to terminate the sandwiching of the preceding CZs (and turn them into CNOTs) as well as restart the sandwiching of the succeeding CZ. So we sandwich all data qubits in this terms (be they L or R) and all check qubits in this check.
  



  if matrix == 'B' or matrix == 'BT':
    # We're doing B or B^T. If check is X i.e. Hx = [A|B] then it's interacting with R data qubits (L-data are idling) and if check is Z i.e. Hz = [B^T|A^T] then it's interacting with L data qubits (R are idling)
    if check == 'X':
      qC = registers.qX
      qD = registers.qR
      qD_idle = registers.qL
      if ONLYCZs == False:
        hadamard(circ, qC + qD, errors)
        idle(circ, qD_idle, idle_during['H'])
        tick(circ)
      elif ONLYCZs == True:
        hadamard(circ, qC + qD + qD_idle, errors) # we will apply the final Hadamards (converting CZs to CNOTs) here on the data qubits that are not being swapped (would usually be idle.)
        tick(circ)
    
    elif check == 'Z':
      qC = registers.qZ
      qD = registers.qL
      qD_idle = registers.qR
      hadamard(circ, qC + qD, errors)
      idle(circ, qD_idle, idle_during['H'])
      tick(circ)


  if matrix == 'A' or matrix == 'AT':
    
    # We're NOW doing A or A^T. 
    # If check is X i.e. Hx = [A|B] then it's interacting with L data qubits (R-data are idling) and if check is Z i.e. Hz = [B^T|A^T] then it's interacting with R data qubits (L are idling)
    if check == 'X':
      qC = registers.qX
      qD = registers.qL
      qD_idle = registers.qR
      if ONLYCZs == False:
        hadamard(circ, qC + qD, errors)
        idle(circ, qD_idle, idle_during['H'])
        tick(circ)
      elif ONLYCZs == True:
        hadamard(circ, qC + qD + qD_idle, errors)
        tick(circ)
    
    elif check == 'Z':
      qC = registers.qZ
      qD = registers.qR
      qD_idle = registers.qL
      hadamard(circ, qC + qD, errors)
      idle(circ, qD_idle, idle_during['H'])
      tick(circ)


def add_hadamards_after_swap_CZ(matrix, thegate, check, registers, circ, errors, idle_during):
  
  if thegate != 'CZ':
    raise ValueError("This is written to go after a CZ gate")
  
  if matrix == 'B' or matrix == 'BT':
    if check == 'X':
      qC = registers.qX
      qD = registers.qR
      qD_idle = registers.qL
    elif check == 'Z':
      qC = registers.qZ
      qD = registers.qL
      qD_idle = registers.qR
    hadamard(circ, qC + qD, errors)
    idle(circ, qD_idle, idle_during['H'])
    tick(circ)

  if matrix == 'A' or matrix == 'AT':

    if check == 'X':
      qC = registers.qX
      qD = registers.qL
      qD_idle = registers.qR

    elif check == 'Z':
      qC = registers.qZ
      qD = registers.qR
      qD_idle = registers.qL

    hadamard(circ, qC + qD, errors)
    idle(circ, qD_idle, idle_during['H'])
    tick(circ)



'''add_swap_lrc
This adds a swap gate after the final 2q interaction for each stabiliser, swapping each check qubit with the last data qubit it interacts with. This means every physical qubit gets measured every other round (while the data is preserved over rounds, unmeasured), providing a mechanism to detect leakage with a leakage-detecting measurement or at least to stop leakage propagating. As per Natalie Brown et al. (10.1088/1367-2630/ab3372) when working with CNOTs the circuit reduces to just adding a reversed-direction CNOT before the final CNOT. 

Optional reading:
  We need to make sure we are inserting a reversed-direction CNOT before what is indeed the final CNOT of each check.
  Sometimes only one of the polynomials have a term for a certain value of j. So . If both A and B have this final j then we will insert a swapped direction CNOT before B's CNOTs. If only one of them then we swap the direction before that one.

  We also need to guarantee that the final data qubits that the Z-checks interact with are of the opposite type to the final data qubits that the X-checks interact with, ensuring that every qubit is measured every other round (rather than accidentally only swapping with one data-qubit type.
  
  In X-checks we apply 2q gates according to polynomial A then polynomial B in Hx = [A|B] for each j value. E.g. if j = 3 and A has x^2 y^3 and B has y^3 we apply the interactions for i = 2 first between X-check qubits and L-data qubits then for i = 0 between X-check qubits and R-data qubits. In Z-checks we apply A^T then B^T in Hz = [B^T | A^T], implying the opposite order: R-data *then* L-data.
  Note we are applying j's in descending order for Z-checks (JTunion in bbparamfuncs.py has been reversed - search JTunion.reverse()), where j_max was the final j-value applied for X-checks, -j_max is the final applied in Z-checks. These both have equivalent values for i that appear in equivalent polynomials. E.g. X-checks might have i_0, i_1 in A and i_2 in B whereas the Z-checks will have -i_0, -i_1 in A^T and -i_2 in B^T. By applying 2q gates in Z-checks in the order A^T then B^T and applying this equivalent j_value last.
  This guarantees that the final data qubits that the Z-checks interact with are of the opposite type to the final data qubits that the X-checks interact with.

  When a specific code is decided upon to be realised in hardware, the j_value order can be chosen in a more optimised fashion to reduce shuttling time, while still ensuring that the final data qubits each check interacts with (and is swapped with) is of an opposite type. For now, where we are simulating multiple codes, we use this method of applying ascending j's for X-checks then descending j's for Z-checks and within each j applying the i values for A then B (or A^T then B^T) implying L then R (or R then L) data qubit interactions.'''
def add_swap_lrc(thegate, jval, theunion, code, circ, registers, errors, idle_during, sequential_gates, check, reuse_check_qubits):


  # a is either A or A^T
  # b is either B or B^T

  if thegate == 'CNOT':
    thegate = 'CX'


  if check == 'X':
    qC = registers.qX
    a_ij = code.Aij
    b_ij = code.Bij
    a = 'A'
    b = 'B'
    as_is = code.As_is
    bs_is = code.Bs_is
  
  elif check == 'Z':
    qC = registers.qZ
    a_ij = code.ATij
    b_ij = code.BTij
    a = 'AT'
    b = 'BT'
    as_is = code.ATs_is
    bs_is = code.BTs_is
  

  in_a_ij = False
  in_b_ij = False
  if jval in as_is:
    in_a_ij = True
  if jval in bs_is:
    in_b_ij = True


  if in_a_ij == True and in_b_ij == True: # i.e. this j value is in both A/A^T and B/B^T

    # Add a's 2q gates:
    add_2q_gates(thegate, a, circ, jval, code, registers, errors, idle_during, sequential_gates) # adds for all i value(s) in A or AT


    # Add b's 2q gates, with a reversed one inserted before the last one
    for idx, ival in enumerate(bs_is[jval]):
      if idx == len(bs_is[jval]) - 1: # we're on the last i value so do the swaperoo

        if thegate == 'CZ':
          add_hadamards_before_swap_CZ(b, thegate, check, registers, circ, errors, idle_during)

        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates, swap_ctrl_target = True)

        if thegate == 'CZ': # The other Hadamard in the sandwich.
          add_hadamards_after_swap_CZ(b, thegate, check, registers, circ, errors, idle_during)

        # Ok, we've now added the reversed CNOT, or the CZs sandwiched by hadamards. 
        # We now add the last timestep's two qubit gates:

        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)

        # and then update the registers because things have swapped around.
        update_qubit_indices(code, registers, b, check, ival, jval, reuse_check_qubits)

      else:
        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)


  if in_a_ij == True and in_b_ij == False:  # This j value only has i value(s) in A or A^T

    for idx, ival in enumerate(as_is[jval]):
      if idx == len(as_is[jval]) - 1: # we're on the last ival so do the swaperoo


        if thegate == 'CZ': # to do a reversed CNOT we sandwich both ends of the CZ by Hadamards. This is on the target in order to make the gate a CNOT. It is also on the control in order to terminate the sandwiching of the preceding CZs (and turn them into CNOTs) as well as restart the sandwiching of the succeeding CZ. So we sandwich all data qubits in this terms (be they L or R) and all check qubits in this check.
          add_hadamards_before_swap_CZ(a, thegate, check, registers, circ, errors, idle_during)

        add_2q_gates_for_this_ij_value(thegate, a, circ, ival, jval, code, registers, errors, idle_during, sequential_gates, swap_ctrl_target = True)
        
        if thegate == 'CZ': # Other hadamard of the sandwich
          add_hadamards_after_swap_CZ(a, thegate, check, registers, circ, errors, idle_during)


        add_2q_gates_for_this_ij_value(thegate, a, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)

        update_qubit_indices(code, registers, a, check, ival, jval, reuse_check_qubits)

      else:
        add_2q_gates_for_this_ij_value(thegate, a, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)

  if in_a_ij == False and in_b_ij == True:  # Only in B or B^T

    for idx, ival in enumerate(bs_is[jval]):
      if idx == len(bs_is[jval]) - 1: # we're on the last ival so do the swaperoo



        if thegate == 'CZ': 
          add_hadamards_before_swap_CZ(b, thegate, check, registers, circ, errors, idle_during)

        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates, swap_ctrl_target = True)
        
        
        if thegate == 'CZ': 
          add_hadamards_after_swap_CZ(b, thegate, check, registers, circ, errors, idle_during)
        
        
        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)

        update_qubit_indices(code, registers, b, check, ival, jval, reuse_check_qubits)

      else:
        add_2q_gates_for_this_ij_value(thegate, b, circ, ival, jval, code, registers, errors, idle_during, sequential_gates)





# Not sure I need to write the below function now ... if the registers are updating then the measures should just update too...
# '''swapped_measure
# The swap leakage reduction circuit of add_swap_lrc has swapped each check qubit with the last data qubit it interacted with. So when measuring between rounds we need to alternate our mesaurements between qubits which were originally data qubits and which were originally check qubits.
# # Note to self, could be good to do this by simply adding an integer to the index of the register. As in, if we were going to measure qZ but now would like to measure qL, then we just write Z + 1. We will also need to track whether the Z-checks or X-checks are swapping with L or R data qubits because in the following rounds we need to reverse the direction of the two-qubit gates.
# '''


'''apply_cyclic_shift_and_2q_gates
Append to a stim circuit for a BB code the required cyclic shifts and two-qubit gates to measure the stabilisers. This is according to Algorithm 2 of Ye ... Delfosse (2 × L array; 2508.01879).
This will also add swap-LRC if global SWAPLRC == True and update the registers qX, qL, qR, qZ depending on which were swapped.
Accepts inputs:
- circ: the stim circuit to be appended to
- jval_prev: the previous arrangement of modules. If jval_prev = j this implies that check qubit module M^a_0 was aligned with data qubit module M^d_((0 + j) % m)
- check: whether we are performing the X-stabiliser or Z-stabiliser checks
- code: class containing paramaters of the BB code, e.g. Hx = [A|B], Hz = [B^T|A^T] etc.
- registers: indices for the stim circuit of check qubits, data qubits
- errors: what errors and probabilities are on each operation
- sequential (bool): whether the two-qubit gates within a module are applied sequentially (in serial) or in parallel
- reuse_check_qubits: whether the one set of check qubits are being reused (possible because we're doing non-interleaved syndrome extraction, i.e. X-checks THEN Z-checks)'''
def apply_cyclic_shift_and_2q_gates(circ, jval_prev, check, code, registers, errors, idle_during, sequential_gates, reuse_check_qubits):

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
        if errors['shift'].p > 0:  
  
          j_dif = abs((jval % m) - (jval_prev % m))


          if j_dif > 0: # i.e. if a shift needs to occur
            
            update_shift_probs(j_dif, errors, idle_during) # updates noise according to length of shift
            
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

        if check == 'X':  # Hx = [A|B]
          
          if ONLYCZs == False:
            thegate = 'CNOT'
          elif ONLYCZs == True:
            thegate = 'CZ'
          
          if SWAPLRC == True: # add_swap_lrc will add a SWAP gate by adding a reversed CNOT before the final CNOT or, with CZs, a CZ with both its target and control sandwiched between hadamards before the final CZ (as well as some added or cancelled hadamards after the final CZ)
            if jval == theunion[-1]: # If we're at the final j value add the reversed CNOTs and the final CNOTs for the final i value.
              add_swap_lrc(thegate, jval, theunion, code, circ, registers, errors, idle_during, sequential_gates, check, reuse_check_qubits)
            else:
              add_2q_gates(thegate, 'A', circ, jval, code, registers, errors, idle_during, sequential_gates)
              add_2q_gates(thegate, 'B', circ, jval, code, registers, errors, idle_during, sequential_gates)

          elif SWAPLRC == False:

            add_2q_gates(thegate, 'A', circ, jval, code, registers, errors, idle_during, sequential_gates)
            add_2q_gates(thegate, 'B', circ, jval, code, registers, errors, idle_during, sequential_gates)



        elif check == 'Z': # Hz = [B^T|A^T]

          if ONLYCNOTs == False:
            thegate = 'CZ'
          elif ONLYCNOTs == True:
            thegate = 'CNOT'
          
          
          if SWAPLRC == True: # Do a reversed CNOT before the final CNOTs in order to implement a swap gate ... 
            if jval == theunion[-1]: #and check == 'X': # TESTING -- STOPPING THIS FROM FIRING # If we're at the final j value add the reversed CNOTs and the final CNOTs for the final i value.
              add_swap_lrc(thegate, jval, theunion, code, circ, registers, errors, idle_during, sequential_gates, check, reuse_check_qubits)
            else:
              add_2q_gates(thegate, 'AT', circ, jval, code, registers, errors, idle_during, sequential_gates)
              add_2q_gates(thegate, 'BT', circ, jval, code, registers, errors, idle_during, sequential_gates)


          elif SWAPLRC == False:

            add_2q_gates(thegate, 'AT', circ, jval, code, registers, errors, idle_during, sequential_gates)
            add_2q_gates(thegate, 'BT', circ, jval, code, registers, errors, idle_during, sequential_gates)


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






'''add_relax_then_leak
   Appends to a deltakit_stim circuit a RELAX gate then a LEAKAGE gate (interpretable by deltakit_stim) to the given register as per the p_relax and p_leak contained in the dictionary errors for the given operation.
   E.g. add_relax_then_leakage('CNOT', circuit, [0, 1, 2], helios_errors(p))

   We append relax THEN leakage. Note that in your errors dictionary, to have P(a qubit leaks) = P(a leaked qubit relaxes) requires p_relax = p_leak / (1 - p_leak) in the error dictionary so when you draw the tree diagram you get P(a qubit leaks) = p_l and P(a qubit relaxes and does not subsequently leak) = p_l'''
def add_relax_then_leak(operation, circ, register, errors: dict):

  p_relax = errors[operation].p_relax
  
  p_leak = errors[operation].p_leak

  if p_relax is not None and p_relax > 0:
    circ.append("RELAX", register, p_relax)
  if p_leak is not None and p_leak > 0:
    #circ.append("LEAKAGE", register, p_leak)
    # temporary for DEP1 leakage:
    circ.append("DEPOLARIZE1", register, 0.75 * p_leak)



''' idle
Adds an idling error to qubits in register. The idling error is of class Error.
E.g. idle(circuit, [0, 1], idle_during['MZ']) '''
def idle(circuit, register, error: Error):

  p = error.p
  
  if p > 0:
    circuit.append(error.op, register, p)

  if LEAKAGE:
    
    p_relax = error.p_relax
    p_leak = error.p_leak

    if p_relax is not None and p_relax > 0:
      circuit.append("RELAX", register, p_relax)
    if p_leak is not None and p_leak > 0:
      # circuit.append("LEAKAGE", register, p_leak)
      # temporary for DEP1 leakage:
      circuit.append("DEPOLARIZE1", register, 0.75 * p_leak)




''' apply_shuttle_error
Applies a depolarising noise channel of strength p to qubits in 'register' and in stim circuit 'circuit'. todo: make more accurate noise model once we have the info'''
def apply_shuttle_error(circuit, register, errors: dict):

    p = errors['shuttle'].p

    if p > 0:

        circuit.append(errors['shuttle'].op, register, p)
      
    if LEAKAGE:
      add_relax_then_leak('shuttle', circuit, register, errors)


''' apply_shift_error
Applies an error to qubits to simulate them undergoing the cyclic shift required to align check and data qubit modules. Contained in errors is the shift constant. The actual value of p can be fed in to represent longer or shorter cyclic shifts'''
def apply_shift_error(circuit, register, errors):
    
    p = errors['shift'].p
    if p > 0:
        circuit.append(errors['shift'].op, register, p) 


''' apply_merge_error
Once check modules have been cyclically shifted to the data qubit module they need to interact with, we simulate merging their coulomb potentials'''
def apply_merge_error(circuit, register, errors: dict):
    p = errors['merge'].p
    if p > 0:
        circuit.append(errors['merge'].op, register, p)

''' apply_split_error
Simulating splitting the coulomb potentials of check qubit and data qubit modules by appling an error probability to them.'''
def apply_split_error(circuit, register, errors: dict):
    p = errors['split'].p
    if p > 0:
        circuit.append(errors['split'].op, register, p)




''' make_loop_body
Constructs the stabiliser extraction round of a memory experiment using a BB code. This round or 'loop_body' can be repeated arbitrarily. Returns the loop_body which is a stim circuit. Note this does not initialise data qubits as it is designed to follow an already-constructed round-0 of stabiliser measurements which is slightly different to the repeated rounds (namely round-0 initialises the data qubits and only has detectors on stabilisers in the same basis as the preserved logical state because the other ones are non-deterministic).
- jval_prev: gives the previous arrangement of modules before starting this repeated / looped section. I.e. if jval_prev = j, this implies check module M^a_w was aligned with data module M^d_((w + j) % m)
- code: paramaters of the BB code. Hx = [A|B], Hz = [B^T,A^T] etc.
- errors: error rates and operations for each circuit operation
- memory_basis: 'Z' implies preserve logical 0, 'X' implies preserve logical plus
- exclude_opposite_basis_detectors: if this is True, when preserving logical 0 (memory_basis = 'Z') then there are no detectors placed on the X-stabiliser measurments (though they are still performed). This is useful if the decoder being used is uncorrelated (i.e. treats X and Z detector graphs separately) as it reduces the size of the detector error model to be fed to it by taking out unused nodes.
- reuse_check_qubits: one register of check qubits of size n/2 -- is possible as we're doing X-checks then Z-checks
- sequential: whether or not the two-qubit gates within a leg are sequential or in parallel'''
def make_loop_body(jval_prev, code, errors, idle_during, registers, memory_basis, reuse_check_qubits, sequential_gates, exclude_opposite_basis_detectors = False):

    n = code.n
    measurements_per_round = 2 * n if LEAKAGE_HERALDS else n

    loop_body = stim.Circuit()

    qX = registers.qX
    qL = registers.qL
    qR = registers.qR
    qZ = registers.qZ  # qZ will equal qX if reuse_check_qubits == True

    idle(loop_body, qL + qR, idle_during['pause']) # This is to simulate pausing rather than applying the syndrome extraction straight away.


    # X-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    qC = qX
   # Initialise check qubits
    init('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['RZ']) # t_init) # idle data qubits
    tick(loop_body)

    # Hadamard check qubits to |+⟩
    hadamard(loop_body, qC, errors)
    if ONLYCZs == False:
      idle(loop_body, qL + qR, idle_during['H']) # t_init) # idle data qubits
    elif ONLYCZs == True: # we are only using CZ gates so Hadamard all data qubits before the X-checks to effectively have CNOTs
      hadamard(loop_body, qL + qR, errors)
    
    tick(loop_body)

    # Do cyclic shifts to required j-valued modules, apply two-qubit gates for stabilisers and return last j position
    jval_prev = apply_cyclic_shift_and_2q_gates(loop_body, jval_prev, 'X', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)
    # Alrighty we've done the X-check CNOTs!

    # Hadamard check qubits (which have already been shuttled back to racetrack in apply_cyclic_shifts_and_stab_interactions)
    hadamard(loop_body, qC, errors)
    if ONLYCZs == False:
      idle(loop_body, qL + qR, idle_during['H']) # t_had) # idle data qubits during hadamard
    elif ONLYCZs == True:
      hadamard(loop_body, qL + qR, errors) # hadamard data qubits as we were using only CZ gates so they need to be hadamarded for X checks.
    tick(loop_body)

    # Measure check qubits
    measure('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['MZ']) # t_meas) # idle data qubits
    

    # We now place detectors on these check qubit measurements. They compare these measurements and the previous round's X-check measurements (even though the first round's measurements might not have detectors on them if we're in memory Z, these are comparing the *measurements* so do not need previous detectors)
    
    if (memory_basis == 'Z' and not exclude_opposite_basis_detectors) or memory_basis == 'X':
            # Append X-check stabiliser detectors:
            for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are the X-check measurements, and compares each of them to the measurement performed n measurements before it if there are no leakage heralds, or 2n before it if there are, as this is the same measurement in the preceding round
                loop_body.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - measurements_per_round)])

    if LEAKAGE_HERALDS: # append leakage heralds to the qubits just measured (qC)
      loop_body.append("HERALD_LEAKAGE_EVENT", qC, 0) 
      for i in reversed(range(1, n//2 + 1)): # append detectors to the heralds
        loop_body.append("DETECTOR", [stim.target_rec(-i)], i)


    if reuse_check_qubits == True: ## If not re-using, then can reset the Z-check qubits in this same time step.
      tick(loop_body)


    # # Z-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    qZ = registers.qZ
    qL = registers.qL
    qR = registers.qR
    qC = registers.qZ 

    # Initialise Z-check qubits
    init('Z', loop_body, qC, errors) # (note qZ = qX if reuse_check_qubits == True)
    idle(loop_body, qL + qR, idle_during['RZ']) # t_init) # idle data qubits
    tick(loop_body)

    # Hadamard check qubits to |+⟩ and IDLE data qubits:
    if ONLYCNOTs == False:
      hadamard(loop_body, qC, errors)
      idle(loop_body, qL + qR, idle_during['H']) # t_init) # idle data qubits
      tick(loop_body)

    # Apply required cyclic shifts and CZ interactions for Z-checks:
    jval_prev = apply_cyclic_shift_and_2q_gates(loop_body, jval_prev, 'Z', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)

    # Now to hadamard the check qubits (they've already been shuttled back into racetrack)
    if ONLYCNOTs == False:
      hadamard(loop_body, qC, errors)
      idle(loop_body, qL + qR, idle_during['H'])  # idle data qubits
      tick(loop_body)

    # Now measure check qubits
    measure('Z', loop_body, qC, errors)
    idle(loop_body, qL + qR, idle_during['MZ']) # t_meas) # idle data qubits
    
    
    # We now place detectors on these check qubit measurements, comparing them and the previous round's Z-check measurements
    
    if (memory_basis == 'X' and not exclude_opposite_basis_detectors) or memory_basis == 'Z':
            # Append Z-check stabiliser detectors:
            for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are now Z-check measurements, and compares each of them to the measurement performed n measurements before it (this is the same measurement in the preceding round)
                loop_body.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - measurements_per_round)])


    if LEAKAGE_HERALDS: # append leakage heralds to the qubits just measured (qC)
      loop_body.append("HERALD_LEAKAGE_EVENT", qC, 0) 
      for i in reversed(range(1, n//2 + 1)): # append detectors to the heralds
        loop_body.append("DETECTOR", [stim.target_rec(-i)], i)

    if reuse_check_qubits == True:
      tick(loop_body)

    return loop_body








''' make_BB_circuit
This function makes a stim circuit realising a memory experiment for any bivariate bicycle (BB) code [2308.07915] according to algorithm 2 of Tham et al.'s "qubit modules" paper [2508.01879]. All the X-checks are performed in parallel, then all the Z-checks. The two-qubit gates proceed according to ascending j, the powers of y in the BB code's matrices A = sum(x^i * y^j), B = sum(x^i * y^j), the two matrices which form the code's parity check matrices Hx = [A|B] and Hz = [B^T|A^T]. All of the required information is contained in the input object 'code'.

Inputs are:
    - code
            An object which contains all the BB code's parameters, such as parity check matrices.
    - p
            Physical error rate
    - errors
            A dictionary that can be made with errors = tham_modules_errors(p) which has an error operation and corresonding probability for each operation in the circuit.
    - idle_during
            A dictionary that can be made with idle_during = tham_modules_idle_errors(p) which has an error operation and corresponding probability for qubits that are idling while other qubits are undergoing each operation in the circuit
    - num_syndrome_extraction_cycles
            How many rounds of stabiliser measurements to perform (including the first round which just encodes the logical state)
    - memory_basis
            Either 'Z' or 'X' to preserve logical 0 or + respectively in all the logical qubits
    - sequential_gates
            Whether disjoint two-qubit gates are sequential (done in successive time steps) or, instead, all in parallel. This affects idling errors applied.
            More explanation: For a given term x^i⋅y^j the most we can possibly do in parallel is all lm X-check qubits to either the lm L or R data qubits and all lm Z-check qubit to the lm opposite type data qubits. So 2lm check qubits connecting to 2lm data qubits. However, we are doing non-interleaved syndrome extraction, that is X-checks *then* Z-checks, so the most that can be done in a single time step is lm check qubits of one type to lm data qubits of one type. (Doing non-interleaved means we can also halve the check-qubit count). Doing them all in a single time step with high fidelity is possible with techniques such as those in [2603.07548]. Alternatively, as in Quantinuum's Helios [2511.05465], gates are done in separate operation zones to maintain high fidelity and reduce crosstalk, meaning only X gates can be done in a single time step, where X is the number of operation zones. This circuit-builder assumes m operation zones, one for each of the m data qubit modules. When setting sequential gates to true, two-qubit gates are done sequentially rather than all in parallel so only m two-qubit gates can be done in each time step.
    - exclude_opposite_basis_detectors
            If this is True, when preserving logical 0 (+) then there are no detectors placed on the X-(Z-) stabiliser measurements (though they are still performed). This is useful if the decoder being used is uncorrelated (i.e. treats X and Z detector graphs separately) as it reduces the size of the detector error model to be fed to it (and thus increases the speed of the simulations) by removing unused detectors.
    - reuse_check_qubits
            If true we have one check register of size l * m rather than two, which is reused for X-checks then Z-checks. This is advised and possible because the algorithm this code implements (algorithm 2 of 2508.01879) does X-checks then Z-checks
    - only_CZs
            If set to true, this creates a circuit that only uses CZs rather than CXs for X-checks and CZs for Z-checks (it does this by Hadamarding all data qubits at the beginning and end of the X-checks (which now feature CZ gates only as well as the errors for cyclically shifting the modules around))
    - only_CNOTs
            If set to true, this creates a circuit using only CNOTs (X-checks have CNOTs from check qubit (in |+⟩ ) to data qubit, Z-checks have CNOTs from data qubit to check qubit (in |0⟩ ). If BOTH only_CNOTs and only_CZs are False then X-checks use CNOTs and Z-checks use CZs.
    - leakage
            Leakage is where the qubit leaves the computational subspace of |0⟩ and |1⟩ and occupies, or is in a superposition with, another energy level, call it |2⟩. Leakage is usually false, meaning no leakage noise will be introduced. Otherwise, when leakage is True, the dictionary of errors (which also contain the leakage rates during different operations) will be used to add leakage noise. For example, if errors = { "H" : Error("DEPOLARIZE1", p_gate, p_leak, p_relax} then after a Hadamard gate and its noise we append RELAX(p_relax) and LEAKAGE(p_leak). These instructions are not recognised by stim but are recognised by deltakit_stim (see examples/leakage_intro_to_deltakit_stim.ipynb) and mean that any previously leaked qubit has a p_relax chance of returning to the computational subspace (where it is a maximally depolarised state) and then it, or any non-leaked qubit, has a p_leak chance of leaking. Leakage is modelled using the depolarising leakage model (10.1088/1367-2630/ab3372) where the qubit maximally depolarises and any qubit it later interacts with also maximally depolarises.
    - leakage_repumping
            As per 10.1103/PhysRevLett.124.170501 (note this is on Ytterbium ions) this is pumps the qubit such that if it is leaked it returns to the qubit manifold as the maximally mixed state with probability 1 - 1/3^n where n is the number of repumping cycles. This imparts a 'memory error' onto non-leaked qubits of n*2*10^-5. For now we will simulate 2 repumping cycles after every two-qubit gate if this is set to true. We also apply the memory error to all qubits whether they are leaked or not.
    - num_repumping_cycles
            As per 10.1103/PhysRevLett.124.170501, given num_repumping_cycles = c, p_relax = 1 - 1/3^c  and p_error = c * 2e-5 . I.e. You reduce the leaked population by a third every time you repump, but introduce an error of 2e-5 on non-leaked qubits (which we apply to all qubits whether they're leaked or not).
    - leakage_heralds 
            If true, the measurements of qubits are also able to return the result 'leaked'. This is simulated by appending the Deltakit Stim gate 'HERALD_LEAKAGE_EVENT()' after each measurement and applying a detector to it.
    - swap-LRC
            "swap leakage reduction circuit". If set to true, an additional CNOT timestep is inserted before the very last CNOT timestep in each of the X-checks and Z-checks. This CNOT interacts the exact same two qubits as the final CNOT it now precedes, but with control and target reversed. The effect of inserting this one additional reversed CNOT is doing the final CNOT and a SWAP gate. So the data qubits and check qubits have swapped roles. Consequently, every qubit will be measured every other round to prevent leakage lasting the whole memory time.
    - check_deteting_regions
            In general this should always be set to True. It uses stim's in-built circuit.detecting_regions() to check if there are any detectors which, in the absence of noise, don't commute and hence produce non-determenistic outcomes. In the absence of noise all detectors should commute so this is a problem if they don't. The only time this should be set to false is when tinkering with the circuit builder and adding and moving detectors, so you might want to look at a visual of a circuit before having finished all the detectors (and consequently want to build the circuit without checking the detecting regions)
    '''
def make_BB_circuit(
    code: Any,
    p: float = 1e-3,
    errors: Any = uniform_errors(1e-3),
    idle_during: Any = uniform_idling(1e-3),
    num_syndrome_extraction_cycles: int = 2,
    memory_basis: str = "Z",
    sequential_gates: bool = False,
    exclude_opposite_basis_detectors: bool = False,
    reuse_check_qubits: bool = False,
    only_CZs: bool = False,
    only_CNOTs: bool = False,
    leakage: bool = False,
    leakage_repumping: bool = False,
    num_repumping_cycles: int = 4,
    leakage_heralds: bool = False,
    swap_LRC: bool = False,
    loss: bool = False,
    check_detecting_regions: bool = True
):
    pass


    if only_CZs and only_CNOTs:
      raise ValueError("Only ONE of only_CZs and only_CNOTs can be True")

    global ONLYCZs # inelegantly using a macro here as it's a late addition (saves me putting it as a new argument in all the sub-functions)
    ONLYCZs = only_CZs
    global ONLYCNOTs
    ONLYCNOTs = only_CNOTs
    global SWAPLRC 
    SWAPLRC = swap_LRC
    global LEAKAGE
    LEAKAGE = leakage
    global LEAKAGE_REPUMPING 
    LEAKAGE_REPUMPING = leakage_repumping
    global LEAKAGE_HERALDS
    LEAKAGE_HERALDS = leakage_heralds
    global REPUMPING_CYCLES
    REPUMPING_CYCLES = num_repumping_cycles


    circ = stim.Circuit()

    registers = make_registers(code.l, code.m, reuse_check_qubits = reuse_check_qubits)
    qX = registers.qX
    qL = registers.qL
    qR = registers.qR
    qZ = registers.qZ  


    add_qubit_coordinates(circ, code, registers, reuse_check_qubits)

   ####### Round 0:  ###############

    # Is different to other rounds in that it:
    # - wait's a time step to initalise data qubits (we are assuming can't prepare |+⟩ state directly so need RZ then H)
    # - in round 0 only puts detectors on X (Z) -checks if preparing logical |+⟩ (|0⟩)


    ############# INITIALISE QUBITS:  !!!!!!!!!!!

    # # # For diagram (get the data qubit resets and measures out of the picture)
    # init('Z', circ, qL + qR, errors)
    # tick(circ)


    # Initialise X-check qubits
    init('Z', circ, qX, errors)

    if memory_basis == 'Z': 
      if ONLYCZs == True:
        init('Z', circ, qL + qR, errors) # need to initialise a step earlier and then hadamard the data qubits  # comment out this line for diagram
    if memory_basis == 'X':
      if ONLYCZs == False:
        init('Z', circ, qL + qR, errors)  # comment out this line for diagram

    tick(circ)


    # Hadamard check qubits to |+⟩
    hadamard(circ, qX, errors)


    # Deal with data qubits:
    if ONLYCZs == True:
      if memory_basis == 'Z':  
        hadamard(circ, qL + qR, errors)

      if memory_basis == 'X': # If only doing CZs, need to Hadamard all the data qubits before the CZs of the X-checks (to make them CNOTs). If memory basis is X though, where you usually prepare in |0⟩ then Hadamard to |+⟩, this means the two hadamards cancel out and all you have to do is prepare in Z here.
        init('Z', circ, qL + qR, errors)  # comment out this line for diagram

    if ONLYCZs == False:
      if memory_basis == 'Z':
        init('Z', circ, qL + qR, errors)   # comment out this line for diagram
      if memory_basis == 'X':
        hadamard(circ, qL + qR, errors)

    tick(circ)


    ######## X-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Junion = code.Junion
    jval_0 = Junion[0] # we assume the starting arrangement of the modules is M^a_w with M^d_((w + j) % m), i.e. no cyclic shift errors initially
    
    # Do cyclic shifts to required j-valued modules (and return last j position)
    # This will also add swap-LRC if swap_LRC == True and update the registers qX, qL, qR, qZ depending on which were swapped.
    jval_prev = apply_cyclic_shift_and_2q_gates(circ, jval_0, 'X', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)



    # Now to hadamard the check qubits (they've already been shuttled back into racetrack in apply_cyclic... function)
    
    if SWAPLRC == False:
      hadamard(circ, qX, errors)
      # Also hadamard data qubits if we were using only CZ gates (sandwiching all the CZ gates with hadamards turns their targets to CNOTs) else idle them ... unless we did a swap_lrc with CZs then it turns out the hadamards cancel
      if ONLYCZs == False:
        idle(circ, qL + qR, idle_during['H']) # t_had)
      elif ONLYCZs == True:
        hadamard(circ, qL + qR, errors)
      tick(circ)

    elif SWAPLRC == True:
      if ONLYCZs == False:
        hadamard(circ, qX, errors)
        idle(circ, qL + qR, idle_during['H'])
        tick(circ)
      elif ONLYCZs == True:
        pass  # hadamards at the end cancel.
    


    # Now measure the check qubits
    
    measure('Z', circ, qX, errors)

    idle(circ, qL + qR, idle_during['MZ']) # t_meas)


    # If preserving logical plus we put detectors on these measurements in the first round:
    n = code.n
    if memory_basis == 'X':
        for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
            circ.append("DETECTOR", [stim.target_rec(-i)])


    if LEAKAGE_HERALDS:  # These are appending correctly, just need to set subsequent detectors to account for if they are here or not 

      circ.append("HERALD_LEAKAGE_EVENT", qX, 0) # set p = 0 i.e. herald tells with certainty if a qubit has leaked rather than potentially making a mistake as to whether it's leaked or not
      
      for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
        circ.append("DETECTOR", [stim.target_rec(-i)], i) # gunna chuck on coordinates ... see what that does



    if reuse_check_qubits == True:
      tick(circ)

    



    # # Z-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    qL = registers.qL  
    qR = registers.qR
    qZ = registers.qZ   ### NEED THIS UPDATE

    # Initialise Z-check qubits
    init('Z', circ, qZ, errors) # (note qZ = qX if reuse_check_qubits == True)
    idle(circ, qL + qR, idle_during['RZ']) # idle data qubits
    tick(circ)

    # Hadamard check qubits to |+⟩ and IDLE data qubits:
    if ONLYCNOTs == False:  # If we're doing Z-checks with CNOTs gates then don't need to Hadamard check qubits, just have reversed CNOTs
      hadamard(circ, qZ, errors)
      idle(circ, qL + qR, idle_during['H']) # idle data qubits
      tick(circ)


    if SWAPLRC:
      old_qZ = qZ.copy() # qZ will be updated in next line so saving this version
      old_qL = qL.copy()
      old_qR = qR.copy()

    # Apply required cyclic shifts and two-qubit interactions for Z-checks:
    jval_prev = apply_cyclic_shift_and_2q_gates(circ, jval_prev, 'Z', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)


    # Now to hadamard the check qubits (they've already been shuttled back into racetrack)
    if ONLYCNOTs == False:
      if SWAPLRC == False:
        hadamard(circ, qZ, errors)
        idle(circ, qL + qR, idle_during['H'])
        tick(circ)
      if SWAPLRC == True: # It turns out that you need to Hadamard what WERE the check qubits during Z-checks with CZ gates
        hadamard(circ, old_qZ, errors)
        idle(circ, old_qL + old_qR, idle_during['H'])
        tick(circ)




    # Now measure check qubits
    measure('Z', circ, qZ, errors)

    idle(circ, qL + qR, idle_during['MZ']) # t_meas)


    # If preserving logical zero we put detectors on these (Z) check measurements in the first round:
    n = code.n
    if memory_basis == 'Z':
        for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
            circ.append("DETECTOR", [stim.target_rec(-i)])


    if LEAKAGE_HERALDS:

      circ.append("HERALD_LEAKAGE_EVENT", qX, 0) # set p = 0 (last argument) i.e. herald tells with certainty if a qubit has leaked rather than potentially making a mistake as to whether it's leaked or not
      
      for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
        circ.append("DETECTOR", [stim.target_rec(-i)], i)


    if reuse_check_qubits == True:
      tick(circ)
    


    # # # for diagram:
    # tick(circ)




    ## SUBSEQUENT ROUNDS:


    ## Make repeated / looped stabiliser measurement rounds:

    if swap_LRC == False:
      if num_syndrome_extraction_cycles > 1:
      
        loop_body = make_loop_body(jval_prev, code, errors, idle_during, registers, memory_basis, reuse_check_qubits, sequential_gates, exclude_opposite_basis_detectors)
    
        # Append loop_body to circuit:
        circ = circ + (num_syndrome_extraction_cycles - 1) * loop_body

    elif swap_LRC == True: # Just manually tack on the extra rounds rather than making a 'repeat' section. Could potentially improve this with a repeated loop_body, however depending on whether we reuse_check_qubits or not it requires different numbers of rounds to be back in original position, so chosen number of stabiliser extraction rounds would need to be the correct number of rounds long. Could write code that just manually tacks on extra rounds, unless num_rounds is a multiple of the correct number, then it could be made into a loop body ... but most of the time apart from specific choices of loops the rounds will just be tacked on anyway. Or there could be a repeated number of rounds in loop body with the extra ones added at the end. Anyway ... having a repeat instruction probably doesn't speed up stim anyway.
      if num_syndrome_extraction_cycles > 1:
        for rep in range(num_syndrome_extraction_cycles - 1):
          
          '''start'''

          offset = n if LEAKAGE_HERALDS else 0

          # Reset X-check qubits
          init('Z', circ, qX, errors)
          idle(circ, qL + qR, idle_during['RZ']) # idle data qubits
          tick(circ)

          # Hadamard X-check qubits to |+⟩ (and Hadamard data qubits if doing only CZ gates)
          hadamard(circ, qX, errors)
          if ONLYCZs == True: # Only doing CZ gates, so also need to Hadamard data qubits (so that the CZs on them will act like CNOTs)
            hadamard(circ, qL + qR, errors)
          elif ONLYCZs == False:
            idle(circ, qL + qR, idle_during['H']) # otherwise just idle data qubits
          tick(circ)

          # Do cyclic shifts and two-qubit gates (accepts previous j_val -- previous position of modules -- and updates it). Might also contain swap-LRC
          jval_prev = apply_cyclic_shift_and_2q_gates(circ, jval_prev, 'X', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)

          # Hadamard check qubits (which have already been shuttled back to racetrack in apply_cyclic_shifts_and_stab_interactions)
          if SWAPLRC == False:
            hadamard(circ, qX, errors)
            if ONLYCZs == False:
              idle(circ, qL + qR, idle_during['H']) # t_had) # idle data qubits during hadamard
            elif ONLYCZs == True:
              hadamard(circ, qL + qR, errors) # hadamard data qubits as we were using only CZ gates so they need to be hadamarded for X checks. 
            tick(circ)
          
          if SWAPLRC == True:
            if ONLYCZs == False:
              hadamard(circ, qX, errors)
              idle(circ, qL + qR, idle_during['H']) # t_had) # idle data qubits during hadamard
              tick(circ)
            elif ONLYCZs == True: # It turns out that when doing the swap LRC with CZ gates on the X checks, the final hadamard of the final CNOT (now a CZ sandwiched by hadamards) cancels with the hadamard just before the check-qubit measurement
              pass


          # Measure check qubits:
          measure('Z', circ, qX, errors)
          idle(circ, qL + qR, idle_during['MZ'])


          # Place detectors on these check qubit measurements which compare them and the previous round's C-check measurements
          n = code.n
          if (memory_basis == 'Z' and not exclude_opposite_basis_detectors) or memory_basis == 'X':
                  # Append X-check stabiliser detectors:
                  for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are the X-check measurements, and compares each of them to the same measurement performed in the round before it: n measurements before it if no leakage heralds, 2n measurements before it if there are (leakage herald adds to measurement record)
                      circ.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - n - offset)]) 

          if LEAKAGE_HERALDS: # append leakage heralds to the X-check qubits just measured
            circ.append("HERALD_LEAKAGE_EVENT", qX, 0) 
            for i in reversed(range(1, n//2 + 1)): # append detectors to the heralds
              circ.append("DETECTOR", [stim.target_rec(-i)], i)


          if reuse_check_qubits: # if not will reset other check qubits in this same time step
            tick(circ)

          ## Z-checks:
          
          # qX = registers.qX; qL = registers.qL; qR = registers.qR; qZ =registers.qZ   # For some reason don't need this update
        
          # Initialise Z-check qubits:
          init('Z', circ, qZ, errors)
          idle(circ, qL + qR, idle_during['RZ'])
          tick(circ)

          # Hadamard check qubits to |+⟩, unless using only CNOT gates then we leave them in |0⟩ and do reversed-direction CNOT gates for the Z parity check
          if ONLYCNOTs == False: 
            hadamard(circ, qZ, errors)
            idle(circ, qL + qR, idle_during['H'])
            tick(circ)

          if SWAPLRC:
            old_qZ = qZ.copy() # qZ will be updated in next line so saving this one 
            old_qL = qL.copy()
            old_qR = qR.copy()


          # Apply required cyclic shifts and CZ interactions for Z-checks, also update registers if required by swaplrc.
          jval_prev = apply_cyclic_shift_and_2q_gates(circ, jval_prev, 'Z', code, registers, errors, idle_during, sequential_gates, reuse_check_qubits)

          # Now to hadamard the check qubits (they've already been shuttled back into racetrack)
          if ONLYCNOTs == False:
            if SWAPLRC == False:
              hadamard(circ, qZ, errors)
              idle(circ, qL + qR, idle_during['H'])
              tick(circ)
            elif SWAPLRC == True: # It turns out that you need to Hadamard what WERE the check qubits
              hadamard(circ, old_qZ, errors)
              idle(circ, old_qL + old_qR, idle_during['H'])
              tick(circ)


          # Now measure check qubits
          measure('Z', circ, qZ, errors)
          idle(circ, qL + qR, idle_during['MZ']) # t_meas) # idle data qubits
          
          
          # We now place detectors on these check qubit measurements, comparing them and the previous round's Z-check measurements
          n = code.n
          if (memory_basis == 'X' and not exclude_opposite_basis_detectors) or memory_basis == 'Z':
                  # Append Z-check stabiliser detectors:
                  for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2]), these are now Z-check measurements, and compares each of them to the measurement performed n measurements before it (this is the same measurement in the preceding round)
                      circ.append("DETECTOR", [stim.target_rec(-i), stim.target_rec(-i - n - offset)])

          if LEAKAGE_HERALDS: # append leakage heralds to the Z-check qubits just measured
            circ.append("HERALD_LEAKAGE_EVENT", qZ, 0) 
            for i in reversed(range(1, n//2 + 1)): # append detectors to the heralds
              circ.append("DETECTOR", [stim.target_rec(-i)], i)


          if reuse_check_qubits == True: # otherwise can do this measurement during same time step as other check qubits are reset
            tick(circ)
          
          '''end'''


    # # Measure all data qubits:
    measure(memory_basis, circ, qL + qR, errors)
    
    # Add leakage heralds onto these measurements:
    if LEAKAGE_HERALDS: # append leakage heralds to the qubits just measured (qL + qR)
      circ.append("HERALD_LEAKAGE_EVENT", qL + qR, 0) 
      for i in reversed(range(1, n + 1)): 
        circ.append("DETECTOR", [stim.target_rec(-i)], i)



    ### Add final detectors:
    add_final_detectors(circ, code, memory_basis)

    ## Annotate logical observables (the Lx's or Lz's if mem X or Z):
    add_logical_observables(circ, code.n, code.Lx, code.Lz, memory_basis)  


    if check_detecting_regions:
      detecting_regions = circ.detecting_regions() # a test to see that it has valid detecting regions (deterministic in the absence of noise)

    # Save circuit:
    # circ.to_file(f"../circuits/nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},b={memory_basis},noise={noise},r={num_syndrome_extraction_cycles},code=BB,l={l},m={m},A='{''.join(str(x) + str(y) for x, y in Aij)}',B='{''.join(str(x) + str(y) for x, y in Bij)}'.stim")

    # # print(circ.to_crumble_url())
    # svg = str(circ.without_noise().diagram("timeline-svg"))
    # with open("output.svg", "w", encoding="utf-8") as f: f.write(svg)

    return circ