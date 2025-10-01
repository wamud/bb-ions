''' circfuncs
Given the parity check matrices (constructed using bbfuncs.py) and paramaters & logical operators (found using bbparams.py) of a Bicycle Bivariate [2308.07915] code, these functions are for constructing a stim circuit that realises a memory experiment using the BB code.
I.e. prepare logical |0⟩ or |+⟩ in all the logical qubits of the BB code, run multiple rounds of stabiliser measurements, measure all the data qubits.
Measurements results from stabilisers and data qubits will then be given to a decoder to see if it correctly predicts what the logical states were.

Note Stim api reference: https://github.com/quantumlib/Stim/wiki/Stim-v1.9-Python-API-Reference'''




''' make_registers
  Makes lists of qubit indices dividing 2n qubits evenly into qA, qB, qC, qD.
  This is most commonly used for a BB code, for example qX, qL, qR, qZ, where qX is the X-check syndrome qubits, qL the 'left' data qubits (appear in left-hand side of Hx; acted on by matrix A), qR the 'right' data qubits (appear in right-hand side of Hx; acted on by matrix B) and qZ, the Z-check syndrome qubits.

  qA: qubits 0 to n/2 - 1
  qB: qubits n/2 to n - 1
  qC: qubits n to 3n/2 - 1 
  qD: qubits 3n/2 to 2n - 1'''
def make_registers(n):
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
  if p != 0:
    circuit.append("Z_ERROR", register, p)

''' initZ
Sets qubits in the list 'register' to the zero state |0⟩.
Adds a reset error (set to |1⟩ ) with probability p'''
def initZ(circuit, register, p = 0):
  circuit.append("R", register)
  if p != 0:
    circuit.append("X_ERROR", register, p)


''' hadamard
Appends a hadamard gate to the stim circuit on qubit(s) specified in 'register'.
After the gate it adds a depolarising noise of strength p (i.e. an error will occur with prob p. Given it occurs, pick one of X, Y or Z at random)'''
def hadamard(circuit, register, p = 0):
  circuit.append("H", register)
  if p != 0:
    circuit.append("DEPOLARIZE1", register, p)


''' idle
Adds uniform depolarising noise of strength p (when an error occurs with probability p pick at random either X, Y or Z)) to qubits in register'''
def idle(circuit, register, p = 0):
  if p != 0:
    circuit.append("DEPOLARIZE1", register, p)

''' idleZ 
Adds dephasing (Z) noise of strength p (a Z operation is applied with probability p) to qubits in 'register. Called "idleZ" as this is used as an idling error but only applies dephasing noise'''
def idleZ(circuit, register, p = 0):
  if p != 0:
    circuit.append("Z_ERROR", register, p)

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