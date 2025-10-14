from bb_ions import *
import stim


# # [[30, 4, 5]] From Ye Delfosse long chain [2503.22071] Table II
# l = 5
# m = 3
# # A = x^0 + x
# # B = x^0 + y + x^2*y^2
# Aij = [(0, 0), (1, 0)]          # the powers (i, j) of x^i * y^j
# Bij = [(0, 0), (0, 1), (2, 2)]
# d = 5

# [[48, 4, 7]] from Ye Delfosse long chain [2503.22071] Table II
#l = 8
#m = 3
## A = x^0 + x
## B = x^0 + y + x^3 * y^2
#Aij = [(0, 0), (1, 0)]
#Bij = [(0, 0), (0, 1), (3, 2)]
#d = 7

# # [[108, 8, 10]] from BB paper [2308.07915] Table III
# l = 9
# m = 6
# # A = x^3 + y + y^2
# # B = y^3 + x + x^2
# Aij = [(3, 0), (0, 1), (0, 2)]
# Bij = [(0, 3), (1, 0), (2, 0)]
# d = 10


# [[144, 12, 12]] 'gross code' from BB paper [2308.07915] Table III
l = 12
m = 6
# A = x^3 + y + y^2
# B = y^3 + x + x^2
Aij = [(3, 0), (0, 1), (0, 2)]
Bij = [(0, 3), (1, 0), (2, 0)]
d = 12


# # [[288, 12, 18]] 'two gross'  from BB paper [2308.07915] Table III
# l = 12
# m = 12
# # A = x^3 + y^2 + y^7
# # B = y^3 + x^1 + x^2
# Aij = [(3, 0), (0, 2), (0, 7)]
# Bij = [(0, 3), (1, 0), (2, 0)]
# d = 18

code = get_code_params(l, m, Aij, Bij, d)



memory_basis = 'Z' # preserve logical 0 or + if basis is Z or X
p = 0.001
noise = 'longchain' # current options are 'longchain' or 'uniform'


num_syndrome_extraction_cycles = code.d_max
reuse_check_qubits = True # one check register, possible as we're doing X-checks then Z-checks
sequential = True # whether or not the two-qubit gates within a leg are sequential or in parallel
exclude_opposite_basis_detectors = True # if this is True, when preserving logical 0 (memory_basis = 'Z') then there are no detectors placed on the X-stabiliser measurments (though they are still performed). This is useful if the decoder being used is uncorrelated (i.e. treats X and Z detector graphs separately) as it reduces the size of the detector error model to be fed to it, taking out unused detectors.


circ = stim.Circuit()


noisetimes = make_longchain_noisetimes(p) if noise == 'longchain' else make_uniform_noisetimes(p) if noise == 'uniform' else make_uniform_agnostic_noisetimes(p) if noise == 'uniform_agnostic' else make_longchain_agnostic_noisetimes(p) if noise == 'longchain_agnostic' else None


registers = make_registers(l, m, reuse_check_qubits = reuse_check_qubits)
qX = registers.qX
qL = registers.qL
qR = registers.qR
qZ = registers.qZ  # qZ will equal qX if reuse_check_qubits == True

t_init = noisetimes.t_init
t_had = noisetimes.t_had
t_merge = noisetimes.t_merge
t_split = noisetimes.t_split
t_cnot = noisetimes.t_cnot
t_cz = noisetimes.t_cz
t_shuttle = noisetimes.t_shuttle
t_shift_const = noisetimes.t_shift_const
t_meas = noisetimes.t_meas
t_idle = noisetimes.t_idle
t_idle_meas = noisetimes.t_idle_meas


## Round 0:
# - wait's a time step to initalise data qubits
# - only puts detectors on X (Z) -checks if preparing logical |+⟩ (|0⟩)


# X-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

qC = qX

# Initialise X-check qubits
init('Z', circ, qC, t_init)
tick(circ)

# Hadamard check qubits to |+⟩ and initalise data qubits to +1 eigenstate of 'memory_basis'
hadamard(circ, qC, t_had)
init(memory_basis, circ, qL + qR, t_init) 
tick(circ)

Junion = code.Junion
jval_0 = Junion[0] # we assume the starting arrangement of the modules is M^a_w with M^d_((w + j) % m), i.e. no cyclic shift errors initially

# Do cyclic shifts to required j-valued modules (and return last j position)
jval_prev = apply_cyclic_shifts_and_stab_interactions(circ, jval_0, 'X', code, registers, noisetimes, sequential)
# Alrighty we've done the X-check CNOTs!


# Now to hadamard the check qubits (they've already been shuttled back into racetrack in apply_cyclic... function)
hadamard(circ, qC, t_had)
idle(circ, qL + qR, t_idle) # t_had)
tick(circ)

# Now measure the check qubits
measure('Z', circ, qC, t_meas)
idle(circ, qL + qR, t_idle_meas) # t_meas)
tick(circ)

# If preserving logical plus we put detectors on these measurements in the first round:
n = code.n
if memory_basis == 'X':
    for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
        circ.append("DETECTOR", [stim.target_rec(-i)])


# # Z-CHECKS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

qC = qZ

# Initialise Z-check qubits
init('Z', circ, qC, t_init) # (note qZ = qX if reuse_check_qubits == True)
idle(circ, qL + qR, t_idle) # t_init) # idle data qubits
tick(circ)

# Hadamard check qubits to |+⟩ and IDLE data qubits:
hadamard(circ, qC, t_had)
idle(circ, qL + qR, t_idle) # t_init) # idle data qubits
tick(circ)

# Apply required cyclic shifts and CZ interactions for Z-checks:
jval_prev = apply_cyclic_shifts_and_stab_interactions(circ, jval_prev, 'Z', code, registers, noisetimes, sequential)

# Now to hadamard the check qubits (they've already been shuttled back into racetrack)
hadamard(circ, qC, t_had)
idle(circ, qL + qR, t_idle) # t_had)
tick(circ)

# Now measure check qubits
measure('Z', circ, qC, t_meas)
idle(circ, qL + qR, t_idle_meas) # t_meas)
tick(circ)

# If preserving logical zero we put detectors on these (Z) check measurements in the first round:
n = code.n
if memory_basis == 'Z':
    for i in reversed(range(1, n//2 + 1)): # appends detectors to last n/2 measurements (i.e. from rec[-1] to rec[-n/2])
        circ.append("DETECTOR", [stim.target_rec(-i)])

## Make repeated / looped stabiliser measurement rounds:
loop_body = make_loop_body(jval_prev, code, noisetimes, registers, memory_basis, reuse_check_qubits, sequential, exclude_opposite_basis_detectors)
# Append to circuit:
circ = circ + (num_syndrome_extraction_cycles - 1) * loop_body

# Final measurement of all data qubits:
measure(memory_basis, circ, qL + qR, t_meas)


### Add final detectors:
add_final_detectors(circ, code, memory_basis)


# Add logical observables (the Lx's or Lz's if mem X or Z):
add_logical_observables(circ, code.n, code.Lx, code.Lz, memory_basis)


detecting_regions = circ.detecting_regions() # a test to see that it has valid detecting regions 

# Save circuit:
circ.to_file(f"../circuits/nkd=[[{code.n}_{code.k}_{code.d_max}]],p={p},b={memory_basis},noise={noise},r={num_syndrome_extraction_cycles},code=BB,l={l},m={m},A='{''.join(str(x) + str(y) for x, y in Aij)}',B='{''.join(str(x) + str(y) for x, y in Bij)}'.stim")

# # print(circ.to_crumble_url())
# svg = str(circ.without_noise().diagram("timeline-svg"))
# with open("output.svg", "w", encoding="utf-8") as f: f.write(svg)
