''' noisefuncs'''

import numpy as np


LEG_SPACING = 5e-3 # [m]

SHUTTLE_DISTANCE = 2.5e-3 # 2.5mm to get in and out of a leg. Calling this a "shuttle" as opposed to the shuttling that occurs during a cyclic shift (which we call a "shift")

SHUTTLE_SPEED = 1 # [m/s]


class Error:
    def __init__(self, operation, p):
        self.op = operation # the error operation
        self.p = p # its probability



''' tham_modules_errors
Defines noise values as per Tham ... Delfosse "qubit modules" [2508.01879] (page 4) which uses Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module (i.e. assuming each module is a long chain) plus a cyclic shift error rate of 30p/100 when shifting modules (to align them).
(Note we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)'''
def tham_modules_errors(p):

    errors = {

        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 10), # note is depolarize (as per longcahin paper) as opposed to X_ERROR
        "RX" : Error("DEPOLARIZE1", p / 10), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", p / 10),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p / 10),
        "MX" : Error("Z_ERROR", p / 10),


        # "Qubit modules" [2508.01879] pg. 4
        "shift" : Error("DEPOLARIZE1", 30 * p / 100),
        "shift_prop_to" : None, # shift_prop_to is used to make errors proportional to the length of the shift. Tham et al. say the noise is independent of the length of the shift so set shift_prop_to to None. (This means circfuncs.update_shift_probs will not change the shift error, making it independent of the length of the shift).


        # Additional for our architecture:
        "shuttle" : Error("DEPOLARIZE1", 0), # we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),

    }
    return errors


''' tham_modules_idle_errors
    Defines noise values as per Tham ... Delfosse "qubit modules" [2508.01879] (page 4) which uses Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module (i.e. assuming each module is a long chain) plus a cyclic shift error rate of 30p/100 when shifting modules (to align them).'''
def tham_modules_idle_errors(p):

    idle_during = {
        
        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 100),
        "RX" : Error("DEPOLARIZE1", p / 100), 
        "H" : Error("DEPOLARIZE1", p / 100),
        "CNOT" : Error("DEPOLARIZE1", p / 100),
        "CZ" : Error("DEPOLARIZE1", p / 100),
        "MZ" : Error("DEPOLARIZE1", 30 * p / 100),
        "MX" : Error("DEPOLARIZE1", 30 * p / 100),

        # "Qubit modules" [2508.01879] pg. 4
        "shift" : Error("DEPOLARIZE1", 30 * p / 100), 
        "shift_prop_to" : None, # i.e. set to None means the the shift idling error will always be the value in the line above rather than this constant multiplied by the length of the shift

        # Additional for our architecture:
        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "pause" : Error("DEPOLARIZE1", 0), # this is an idling error applied at the beginning of each round of stabiliser measurements; simulates waiting before each round of stab. measurement

    }

    return idle_during


''' helios_errors
Defines noise values as per Quantinuum's Helios quantum computer [2511.05465]. For a breakdown of how these have been calculated see https://docs.google.com/spreadsheets/d/1WdbadMM03gGbK52di-t6eae_xXPnevyxqVAySLAWTwM/edit?usp=sharing
Note that they combined all transport (cylic shift, merge/split Coulomb potentials, junction enter/exit) and cooling into one "depth-n memory error" (we take depth-1 even though this is an overestimate for our more organised circuit as compared to their random pairings of 98 qubits), which we will divide into two independent shuttle errors (setting P(exactly one Z error from two shuttles) = P(one Z error from depth-1 memory time)) because this is equivalent for all the cyclic shifts while also makes a roughly half as likely Z error for when just doing the final shuttle of qubits before measurements.

Note also that setting p = 0.001 will give all the values as per Helios, most notably a two-qubit gate error rate of 7 × 10^-4 ≈ 1 × 10^-3 and a reset / measure error rate of 1e-3'''
def helios_errors(p):

    helioserrors = {

        # These values are equal to Helios values when p = 0.001:

        "RZ" : Error("X_ERROR", p), 
        "RX" : Error("Z_ERROR", p),  
        "H" : Error("DEPOLARIZE1", 1.4e-2  * p), # Helios value is 1.4e-5, so add a power of 3 to it to account for p = 1e-3
        "CNOT" : Error("DEPOLARIZE2", 7e-1 * p), # 
        "CZ" : Error("DEPOLARIZE2", 7e-1 * p),
        "MZ" : Error("X_ERROR", p),
        "MX" : Error("Z_ERROR", p),

        # Additional for our architecture (all accounted for in shuttle error)
        "shuttle" : Error("DEPOLARIZE1", 12 * p / 100), # we usually define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules, as distinct from the cyclic shift of modules around the racetrack). For Helios noise though it makes more sense to just put other transport errors to zero and just make two shuttles represent the split, shuttle, cyclic shift, shuttle, merge and cooling. That's because usually the process goes
        # Shuttle qubits into leg, merge their coulomb potentials, perform required two qubit gates (all powers of i for that power of j in the BB code's polynomial Σ_{i,j}(x^iy^j) ), split their coulomb potentials, shuttle, cyclic shift to next power of j, repeat. 
        # However the error from Helios incoporates all of that plus cooling! (In fact it incorporates more than all of that as it is the transport error from randomly pairing 98 qubits and perforing two qubit gates on them, whereas our circuit is simpler as the qubits do not require as much rearrangement, so in reality this is an overestimate).
        # So to incorporate that we'll just make two shuttle errors equivalent to all their "depth-1" transport (so as to let the code do its thing at the very beginning and end, saying there is half a depth-1 transport when just shuttling out at the very end before measurement)
        
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, # shift_prop_to is used to make errors proportional to the length of the shift. We are OVERESTIMATING the shift error by using Helios' combined value for shift, merge/split/ junction enter exit and cooling operations from a 98 qubit program that had to randomly sort all the qubits because our program is more organised and just shifts a module of qubits to another module. Interestingly, this Helios overestimate comes out as only slightly less than the Tham et al. estimate of 30p/100
        
    }
    return helioserrors


''' helios_idle_errors
    Defines noise values as per Quantinuum's Helios quantum computer [2511.05465].
    This calculation uses that the idling error is depolarising noise. It is also identical to the cyclic shift error so probably an overestimate as these idling errors will sometimes apply to qubits that aren't moving (however this opens up the simulation possibility of equivalently just saying they're being shuttled around and cooled whenever they're idling).
    For more details on the calculation see https://docs.google.com/spreadsheets/d/1WdbadMM03gGbK52di-t6eae_xXPnevyxqVAySLAWTwM/edit?usp=sharing
    but basically they did an experiment where they randomly pair 98 qubits which are found in a loop. They then have to cyclically shift the loop back and forth until they have the correct pairings in the legs (which can only fit eight qubits at a time and there are two legs) and then perform four two-qubit gates, then intrazone shift, then the other four two-qubit gates after some cooling. They did all the required Coulomb potential split/merge, shifts, shuttles, cooling, junction exit / enter etc. as if they were performing these two-qubit gates but without actually performing the gates and found that the average infidelity on the qubits was 1.6 × 10^-4. They say this is due to spatiotemporal inhomogeneities in the magnetic field, implying it is dephasing noise as this changing B field would slightly change the energy levels of the hyperfine qubit that they're using and thus cause a mismatch between their phase tracking and the natural time evolution of the state under the schrodinger equation (which causes it to precess around the Bloch sphere with a relative phase e^-it(E_0 - E_1)/ℏ). With a calculation as shown in the link we can get p_z of 2.4 × 10^-4 from this average infidelity. We can then reverse calculate an effective T_2 (making the assumption that even un-moving qubits suffer from the same level of dephasing -- an overestimate) using the p_z from a dephasing channel p_z = (1/2)(1 - e^(-t/T_2)) and get T2 = 115s.'''
def helios_idle_errors(T2 = 115):

    # t_2q = 70e-6  # from paper. It's ≈ 70μs
    t_2q = 650e-6  # Have added idling time for a four-ion shift of 280μs and 300μs of cooling and the 70μs two-qubit gate, i.e. 650μs. This simulates being able to do four-ion shifts between two-qubit gates (i.e. simulating just one operation zone when sequential_gates = True in the circuit! (Also it's still an order of magnitude lower than tham_modules idling errors)).
    t_1q = 70e-6    # assuing same as two-qubit gate (overestimate. Usually is at least an order of magnitude faster)
    
    
    ## Actually doesn't matter about measure and reset times for idling as they're dominated by crosstalk errors (which are over an order of magnitude larger)
    # t_m  = 240e-6    # double the H1 quantinuum quantum computer (https://doi.org/10.1038/s41586-021-03318-4) in their extended data Fig. 1. c) (there was no reported time in Helios)
    # t_r  = 310e-6   # reset is measurement + single qubit X
    
    
    helios_idle_during = {

        # Operation : What qubits suffer that are NOT undergoing the operation (i.e. idling while other qubits have that operation done)
        
        "H" : Error("Z_ERROR", p_idle_dephasing(t_1q, T2)),    
        "CNOT" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)), 
        "CZ" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)),   
        
        # Crosstalk errors
        "RZ" : Error("DEPOLARIZE1", 6e-5),
        "RX" : Error("DEPOLARIZE1", 6e-5), 
        "MZ" : Error("DEPOLARIZE1", 6e-5),
        "MX" : Error("DEPOLARIZE1", 6e-5),

        "shuttle" : Error("Z_ERROR", 1.2e-4), 


        # All the below are accounted for in shuttle
        "merge" : Error("Z_ERROR", 0),
        "split" : Error("Z_ERROR", 0),
        "shift" : Error("Z_ERROR", 0), 
        "shift_prop_to" : None,

        "pause" : Error("Z_ERROR", 0),
    }

    return helios_idle_during



''' zero_errors
    Set all errors to zero'''
def zero_errors():
    """
    Defines zero error rates.
    """
    errors = {

        "RZ" : Error("DEPOLARIZE1", 0), # note is depolarize (as per longcahin paper) as opposed to X_ERROR

        "RX" : Error("DEPOLARIZE1", 0), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", 0),
        "CNOT" : Error("DEPOLARIZE2", 0),
        "CZ" : Error("DEPOLARIZE2", 0),
        "MZ" : Error("X_ERROR", 0),
        "MX" : Error("Z_ERROR", 0),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_prop_to" : 0, 
        
    }
    
    return errors



def zero_idling():
    """
    Defines default idling gates and error rates, i.e. the operation and probability of that operation applied to qubits that are idling in a timestep that other qubits are experiencing the key operation.
    """
    idle_during = {
        "RZ" : Error("DEPOLARIZE1", 0),
        "RX" : Error("DEPOLARIZE1", 0), 
        "H" : Error("DEPOLARIZE1", 0),
        "CNOT" : Error("DEPOLARIZE1", 0),
        "CZ" : Error("DEPOLARIZE1", 0),
        "MZ" : Error("DEPOLARIZE1", 0),
        "MX" : Error("DEPOLARIZE1", 0),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),
        "pause" : Error("DEPOLARIZE1", 0),

        "shift_prop_to" : 0,
    }

    return idle_during


''' uniform_errors
Defines a standard depolarising noise channel, as used in original Bravyi et al. BB paper [2308.0791] (see page 16).'''
def uniform_errors(p):

    errors = {

        "RZ" : Error("X_ERROR", p), 
        "RX" : Error("Z_ERROR", p),
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p),
        "MX" : Error("Z_ERROR", p),


        # Qubit module errors -- None
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, 

        # Additional for our architecture - None
        "shuttle" : Error("DEPOLARIZE1", 0), 
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
    }
    
    return errors


def uniform_idling(p):

    idle_during = {
        "RZ" : Error("DEPOLARIZE1", p),
        "RX" : Error("DEPOLARIZE1", p), 
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE1", p),
        "CZ" : Error("DEPOLARIZE1", p),
        "MZ" : Error("DEPOLARIZE1", p),
        "MX" : Error("DEPOLARIZE1", p),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_prop_to" : 0,

        "pause" : Error("DEPOLARIZE1", 0),

    }

    return idle_during


''' p_idle_dephasing
If a qubit is idling for time t, this function returns what the probability of a Z-error occurring on it will be. It takes as inputs t, the time the qubit is idling for, and T2, the characteristic time for dephasing of an idling qubit. For example p_idle(t = 100e-6, T2 = 50s) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience a Z-error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2))'''
def p_idle_dephasing(t, T2):
    p = 0.5 * (1 - np.exp(-t / T2))
    return p

''' dephasing_idle_errors
Defines idle errors as the realistic dephasing noise for a given T2 value with times set by the sources listed below in comments'''
def dephasing_idle_errors(T2):

    t_1q = 7.5e-6     # [27] in Bruzewicz et al. Table 1
    t_2q = 100e-6     # [27] in Bruzewicz et al. Table 1
    t_m  = 200e-6     # Myerson et al. https://doi.org/10.1103/PhysRevLett.100.200502
    t_r  = 207.5e-6   # reset is measurement + single qubit X
    
    idle_during = {
        
        "RZ" : Error("Z_ERROR", p_idle_dephasing(t_r, T2)),
        "RX" : Error("Z_ERROR", p_idle_dephasing(t_r, T2)), 
        "H" : Error("Z_ERROR", p_idle_dephasing(t_1q, T2)),    
        "CNOT" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)), 
        "CZ" : Error("Z_ERROR", p_idle_dephasing(t_2q, T2)),   
        "MZ" : Error("Z_ERROR", p_idle_dephasing(t_m, T2)),
        "MX" : Error("Z_ERROR", p_idle_dephasing(t_m, T2)),

        "shuttle" : Error("Z_ERROR", p_idle_dephasing(SHUTTLE_DISTANCE / SHUTTLE_SPEED, T2)), 
        
        "merge" : Error("Z_ERROR", 0),
        "split" : Error("Z_ERROR", 0),
        
        "shift" : Error("Z_ERROR", 0.1), # will be updated by circfuncs.update_shift_prob as long as p != 0.

        "shift_prop_to" : T2, # T2 time

        "pause" : Error("Z_ERROR", 0),
    }

    return idle_during


def our_uniform_plus_shift_and_shuttle(p, T2):

    SHUTTLE_SPEED = 1 # [m/s]

    errors = {
        
        # Uniform errors:
        "RZ" : Error("X_ERROR", p), 
        "RX" : Error("Z_ERROR", p),
        "H" : Error("DEPOLARIZE1", p),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p),
        "MX" : Error("Z_ERROR", p),

        # Qubit module errors
        "shift_prop_to" : T2, 
        "shift" : Error("Z_ERROR", 0.1), # will be updated each shift
 

        # Additional for our architecture 
        "shuttle" : Error("Z_ERROR", p_idle_dephasing(SHUTTLE_DISTANCE / SHUTTLE_SPEED, T2)), 

        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
    }
    
    return errors






''' idle
Adds an idling error to qubits in register. The idling error is of class Error.
E.g. idle(circuit, [0, 1], idle_during['MZ']) '''
def idle(circuit, register, error: Error):

  p = error.p
  
  if p > 0:
    circuit.append(error.op, register, p)


''' apply_shuttle_error
Applies a depolarising noise channel of strength p to qubits in 'register' and in stim circuit 'circuit'. todo: make more accurate noise model once we have the info'''
def apply_shuttle_error(circuit, register, errors: dict):

    p = errors['shuttle'].p

    if p > 0:

        circuit.append(errors['shuttle'].op, register, p)



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
