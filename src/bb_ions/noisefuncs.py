''' noisefuncs'''

import numpy as np


LEG_SPACING = 5e-3 # [m]

SHUTTLE_DISTANCE = 2.5e-3 # 2.5mm to get in and out of a leg. Calling this a "shuttle" as opposed to the shuttling that occurs during a cyclic shift (which we call a "shift")

SHUTTLE_SPEED = 1 # [m/s]


class Error:
    def __init__(self, operation, p, p_leak = None, p_relax = None):
        self.op = operation # the error operation
        self.p = p # its probability
        self.p_leak = p_leak # the probability of a qubit leaking during this operations
        self.p_relax = p_relax # the propability of a previously-leaked qubit relaxing during this operation



''' longchain_errors
Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module 
(Note we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)'''
def longchain_errors(p):

    errors = {

        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 10), # note is depolarize (as per longcahin paper) as opposed to X_ERROR
        "RX" : Error("DEPOLARIZE1", p / 10), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", p / 10),
        "CNOT" : Error("DEPOLARIZE2", p),
        "CZ" : Error("DEPOLARIZE2", p),
        "MZ" : Error("X_ERROR", p / 10),
        "MX" : Error("Z_ERROR", p / 10),


        # NO shuttling errors (everything below is zero)

        
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, # shift_prop_to is used to make errors proportional to the length of the shift. Tham et al. say the noise is independent of the length of the shift so set shift_prop_to to None. (This means circfuncs.update_shift_probs will not change the shift error, making it independent of the length of the shift).
        "shuttle" : Error("DEPOLARIZE1", 0), # we define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules. For the "shuttling" required for the cyclic shifts we call this "shift" error)
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),


        # Added repumping error from Honeywell:
        "error_per_repump" : Error("DEPOLARIZE1", 2e-2 * p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 

    }
    return errors


''' longchain_idle_errors
    Ye & Delfossse "long chains of trapped ions" [2503.22071] noise values for within each module (i.e. assuming each module is a long chain) '''
def longchain_idle_errors(p):

    idle_during = {
        
        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 100),
        "RX" : Error("DEPOLARIZE1", p / 100), 
        "H" : Error("DEPOLARIZE1", p / 100),
        "CNOT" : Error("DEPOLARIZE1", p / 100),
        "CZ" : Error("DEPOLARIZE1", p / 100),
        "MZ" : Error("DEPOLARIZE1", 30 * p / 100),
        "MX" : Error("DEPOLARIZE1", 30 * p / 100),

        # NO shuttling errors: 
        "shift" : Error("DEPOLARIZE1", 0), 
        "shift_prop_to" : None, # i.e. set to None means the the shift idling error will always be the value in the line above rather than this constant multiplied by the length of the shift
        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "pause" : Error("DEPOLARIZE1", 0), # this is an idling error applied at the beginning of each round of stabiliser measurements; simulates waiting before each round of stab. measurement
    }

    return idle_during



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

        # Added repumping error from Honeywell in case it gets called:
        "error_per_repump" : Error("DEPOLARIZE1", 2e-2 * p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 

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


''' walking_cat_approx
    IonQ's 'walking cat architecture' simulation noise model https://arxiv.org/abs/2604.19481v1
    " the (0, i, j) cyclic gate, which requires im + j transport steps". They also assume fast highway travel for long-ring shift (swapping L and R data qubits being aligned with X or Z checks), e.g. (1, i, j). Each transport step takes 1/20 POC so has p = p/2000. An approximation of the idling error for transport for the gross code is presented below.'''
def walking_cat_errors(p, p_leak = 0):


    # For the gross code we have l = 12, m = 6 and x^3, y, y^2 and y^3, x, x^2 so steps of (Fig. 13 in walking cat paper helpful with this)
    # 3*m = 18, then back 18 plus 1, then another 1 for X-checks. Then 1 step to swap to R data qubits, then 3, then 3 back and 6 across, then another 6 across. Similar for Z checks. So if we just overestimate and let's say 20 transport steps between each power of j. That's 1 POC. Let's overestimate and say each cyclic shift takes two shuttles and each is 20 transport steps so p / 100 -- i.e. normal idling.

    p_relax = round_sig_fig(p_leak / (1 - p_leak), 6)

    errors = {

        "CNOT" : Error("DEPOLARIZE2", p, p_leak, p_relax),
        "CZ" : Error("DEPOLARIZE2", p, p_leak, p_relax),

        "RZ" : Error("DEPOLARIZE1", p / 10, p_leak, p_relax), # note is depolarize (as per longcahin paper) as opposed to X_ERROR
        "RX" : Error("DEPOLARIZE1", p / 10, p_leak, p_relax), # as opposed to Z_ERROR
        "H" : Error("DEPOLARIZE1", p / 10, p_leak, p_relax),
        "MZ" : Error("X_ERROR", p / 10, p_leak, p_relax),
        "MX" : Error("Z_ERROR", p / 10, p_leak, p_relax),

        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, 

        "shuttle" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),

        # Added repumping error from Honeywell in case it gets called:
        "error_per_repump" : Error("DEPOLARIZE1", 2e-2 * p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 

    }
    return errors


''' walking_cat_idle_approx
    IonQ's 'walking cat architecture' simulation noise model https://arxiv.org/abs/2604.19481v1
    Cyclic shift by x^i y^j takes i*m + j time steps. If swapping L and R data qubits it's an extra 1 times step.'''
def walking_cat_idle_approx(p, p_leak = 0):
    
    
    p_relax = round_sig_fig(p_leak / (1 - p_leak), 6)
    
    
    idle_during = {
        
        # Longchain [2503.22071] operations:
        "RZ" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
        "RX" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax), 
        "H" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
        "CNOT" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
        "CZ" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
        "MZ" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),  # different to Tham modules (which had 30p/100 for longer measurement)
        "MX" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),  #                              ditto

        
        "shift" : Error("DEPOLARIZE1", 0), 
        "shift_prop_to" : None, 

        # Additional for our architecture:
        "shuttle" : Error("DEPOLARIZE1", p / 100, p_leak, p_relax),
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

    # Below is distribution of Pauli errors for Rzz gate taken from Table A4 of Helios [2511.05465]: (things that can be conjugated by ZZ to eachother are symmetrical - same orbit of R_ZZ -- as are things that you can swap X and Y to get to eachother (search "gauge freedom" in the paper, this is because the two qubit gate is RZZ(π/2)). I would assume that these errors include the errors introduced by the additional 300

    # IX, IY, ZX, ZY = 4.5e-5
    # XI, YI, XZ, YX = 6e-5
    # XX, YY, XY, YX = 0 # is two orders of magnitude less so might as well be 0
    # IZ, ZI = 1.9e-4
    # ZZ = 6e-5

    rzzprobs = [
        4.5e-5,   # IX
        4.5e-5,   # IY
        19e-5,    # IZ
        5.8e-5,   # XI
        0,        # XX
        0,        # XY
        5.8e-5,   # XZ
        5.8e-5,   # YI
        0,        # YX
        0,        # YY
        5.8e-5,   # YZ
        19e-5,    # ZI
        4.5e-5,   # ZX
        4.5e-5,   # ZY
        6e-5      # ZZ
    ]

    # Recall: Error(error_operation, p_error, p_leak, p_relax)

    helios_errors = {


        # These values are equal to Helios values when input p is 0.001:

        "RZ" : Error("X_ERROR", p, 0, 0), # p_leak = p_relax = 0 because "Typically, ions are initialized using optical pumping techniques which do not result in leakage (https://doi.org/10.1103/PhysRevA.100.032325)"
        "RX" : Error("Z_ERROR", p, 0, 0),  

        # To make these values equal the Helios values (we have constants multiplied by 10^3 so that when input p is 1e-3 they equal the Helios values)
        
        
        "H" : Error("DEPOLARIZE1", 1.4e-2  * p,  1.1e-2 * p,  round_sig_fig( (1.1e-2 * p / (1 - 1.1e-2 * p)), 5)), 
        
        
        # Make p_relax = p_leak / (1 - p_leak) so when you apply relax(p_relax) then leakage(p_leak) and find the joint probabilities, you actually have P(leaked qubit relaxes) = p_leak and P(a qubit leaks) = p_leak

        # For 2q gates leakage from 2QCB (what we use to get the partial pauli error model) is 1.14×10^(−4) (page 8 of Helios paper v1)
        # ⇒ P(at least one leakage) =  1.14 × 10^(−4) = p_l^2+2p_n p_l = p_l^2 + 2(1 - p_l)p_l
        # or 1 − p_n^2 = 1.14 × 10^−4 
        # ⇒ p_l = 5.7 × 10^(−5)  where p_l is the probability of a single qubit leaking so will be applied to each qubit in the gate.
        "CNOT" : Error("DEPOLARIZE2", 7e-1 * p, 5.7e-2 * p, round_sig_fig(5.7e-2 * p / (1 - 5.7e-2 * p), 5)), 


        "CZ" : Error("PAULI_CHANNEL_2", [1e3 * p * prob for prob in rzzprobs], 5.7e-2 * p, round_sig_fig(5.7e-2 * p / (1 - 5.7e-2 * p), 5)), 
        
        
        "MZ" : Error("X_ERROR", p, round_sig_fig(4.2 * p, 5), round_sig_fig(4.2 * p / (1 - 4.2 * p), 5)), 
        "MX" : Error("Z_ERROR", p, round_sig_fig(4.2 * p, 5), round_sig_fig(4.2 * p / (1 - 4.2 * p), 5)), 

        
        # Additional for our architecture (all accounted for in shuttle error)
        

        "shuttle" : Error("Z_ERROR", 1.2e-1 * p, 2.2e-1 * p, round_sig_fig(2.2e-1 * p / (1 - 2.2e-1 * p), 5)), 
        # p_leak = 4.4e-4 (Table A5 Helios paper) during one depth-1 transport. As explained in next comment we're dividing this into two "shuttles" so 2.2e-4 each (which will be the value when p = 1e-3).
        # We usually define "shuttling" as the steps aligning modules before or after they have been cyclically shifted (getting them from the racetrack loop of check qubit modules into the legs that contain the data qubit modules, as distinct from the cyclic shift of modules around the racetrack). For Helios noise though it makes more sense to just put other transport errors to zero and just make two shuttles represent the split, shuttle, cyclic shift, shuttle, merge and cooling. That's because usually the process goes
        # Shuttle qubits into leg, merge their coulomb potentials, perform required two qubit gates (all powers of i for that power of j in the BB code's polynomial Σ_{i,j}(x^iy^j) ), split their coulomb potentials, shuttle, cyclic shift to next power of j, repeat. 
        # However the error from Helios incoporates all of that plus cooling! (In fact it incorporates more than all of that as it is the transport error from randomly pairing 98 qubits and perforing two qubit gates on them, whereas our circuit is simpler as the qubits do not require as much rearrangement, so in reality this is an overestimate).
        # So to incorporate that we'll just make two shuttle errors equivalent to all their "depth-1" transport (so as to let the code do its thing at the very beginning and end, saying there is half a depth-1 transport when just shuttling out at the very end before measurement)
        
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, # shift_prop_to is used to make errors proportional to the length of the shift. We are OVERESTIMATING the shift error by using Helios' combined value for shift, merge/split/ junction enter exit and cooling operations from a 98 qubit program that had to randomly sort all the qubits because our program is more organised and just shifts a module of qubits to another module. Interestingly, this Helios overestimate comes out as only slightly less than the Tham et al. estimate of 30p/100
        
        "error_per_repump" : Error("DEPOLARIZE1", 2e-2 * p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 
    }
    return helios_errors


''' helios_idle_errors
    Defines noise values as per Quantinuum's Helios quantum computer [2511.05465]. When input p = 1e-3 it is exactly as per Helios. 2e-3 is double etc.
    The idling calculation uses that the idling error is dephasing noise and is based on when qubits are idling during their cyclic shift.It is consequently likely an overestimate as these idling errors will sometimes apply to qubits that aren't moving (however this opens up the simulation possibility of equivalently just saying they're being shuttled around and cooled whenever they're idling).
    For more details on the calculation see https://docs.google.com/spreadsheets/d/1WdbadMM03gGbK52di-t6eae_xXPnevyxqVAySLAWTwM/edit?usp=sharing
    Basically, Quantinuum did an experiment where they randomly pair 98 qubits which are found in a loop. They then have to cyclically shift the loop back and forth until they have the correct pairings in the legs (which can only fit eight qubits at a time and there are two legs) and then perform four two-qubit gates, then intrazone shift, then the other four two-qubit gates after some cooling. They did all the required Coulomb potential split/merge, shifts, shuttles, cooling, junction exit / enter etc. as if they were performing these two-qubit gates but without actually performing the gates and found that the average infidelity on the qubits was 1.6 × 10^-4. They say this is due to spatiotemporal inhomogeneities in the magnetic field, implying it is dephasing noise as this changing B field would slightly change the energy levels of the hyperfine qubit that they're using and thus cause a mismatch between their phase tracking and the natural time evolution of the state under the schrodinger equation (which causes the qubit state to precess around the Bloch sphere with a relative phase e^-it(E_0 - E_1)/ℏ). With a calculation as shown in the link we can get p_z of 2.4 × 10^-4 from this average infidelity. We can then reverse calculate an effective T_2 (making the assumption that even un-moving qubits suffer from the same level of dephasing -- an overestimate) using the p_z from a dephasing channel p_z = (1/2)(1 - e^(-t/T_2)) and get T2 = 115s.
    
    LEAKAGE:
    Fig. 5 b) (1QRB), Fig. 6 b) (2QRB) and Fig. 8 b) and c) (Depth-1 Transport) of the Helios paper all show a linear leakage function in sequence length / transport depth so I'll do a linear function in time for leakage on idling qubits. We will assume a non-moving idling qubit suffers the same leakage as the idling transported qubits.
    Calculation:

    - Leakage after 55ms is 4.4e-4  (Table A5)  (Depth-1 Transport)
    ⇒  p_l =  4.4e-4  when t = 0.055s 
    ⇒  m   =  4.4e-4 / 0.055 = 8e-3   
	⇒  p_l =  8e-3 * t                                       <----  idling leakage function'''
def helios_idle_errors(p):


    multiple = p/(1e-3) # When input p = 1e-3 you get exactly the Helios errors. p = 2e-3 would be double etc.
    T2 = 115 # see note above

    m = 8e-3 # from p_l = mt -- the linear leakage function in comments above

    # t_2q = 70e-6  # from paper. It's ≈ 70μs
    t_2q = 650e-6  # Have added idling time for a four-ion shift of 280μs and the 300μs of cooling after a four-ion shift and the 70μs of the actual two-qubit gate, i.e. 650μs. This simulates being able to do four-ion shifts between two-qubit gates (i.e. simulating just one operation zone when sequential_gates = True in the circuit).
    t_1q = 70e-6    # assuming same as two-qubit gate (overestimate so will produce larger idling errors. Usually is at least an order of magnitude faster)
    
    
    
    t_m  = 240e-6    # double the H1 quantinuum quantum computer (https://doi.org/10.1038/s41586-021-03318-4) in their extended data Fig. 1. c) (there was no reported time in Helios)
    t_r  = 310e-6   # reset is measurement + single qubit X

    t_4s = 280e-6 # Time of a four-ion shift
    
    
    helios_idle_during = {

        # Operation : What qubits suffer that are NOT undergoing the operation (i.e. idling while other qubits have that operation done)

        "H" : Error("Z_ERROR", multiple * p_idle_dephasing(t_1q, T2), multiple * m * t_1q, round_sig_fig( multiple * m * t_1q / (1 - multiple * m * t_1q), 5)),  

        "CNOT" : Error("Z_ERROR", multiple * p_idle_dephasing(t_2q, T2), multiple * m * t_2q, round_sig_fig( multiple * m * t_2q / (1 - multiple * m * t_2q), 5)),
        "CZ" :   Error("Z_ERROR", multiple * p_idle_dephasing(t_2q, T2), multiple * m * t_2q, round_sig_fig( multiple * m * t_2q / (1 - multiple * m * t_2q), 5)),   
        
        # Crosstalk errors and leakage idling ( note t_m and t_r not used for the idling on other qubits during reset and measure as crosstalk dominates (over an order of magnitude larger) but are used for the linear leakage function rate)

        # average (including on worst-affected qubits, though we don't have any of those as neighbouring qubits will always also be being measured) crosstalk during MCMR in Helios is 6e-5. Divide this between M and R equally gives p_error = 3.00005e-5
        "MZ" : Error("DEPOLARIZE1", multiple * 3e-5, multiple * m * t_m, round_sig_fig( multiple * m * t_m / (1 - multiple * m * t_m), 5) ), 
        "MX" : Error("DEPOLARIZE1", multiple * 3e-5, multiple * m * t_m, round_sig_fig( multiple * m * t_m / (1 - multiple * m * t_m), 5) ), 
        
        "RZ" : Error("DEPOLARIZE1", multiple * 3e-5, multiple * m * t_r, round_sig_fig( multiple * m * t_r / (1 - multiple * m * t_r), 5) ),
        "RX" : Error("DEPOLARIZE1", multiple * 3e-5, multiple * m * t_r, round_sig_fig( multiple * m * t_r / (1 - multiple * m * t_r), 5) ),

        

        "shuttle" : Error("Z_ERROR", multiple * 1.2e-4, multiple * 2.2e-4, round_sig_fig( multiple * 2.2e-4 / (1 - multiple * 2.2e-4), 5) ), 


        # All the below are accounted for in shuttle
        "merge" : Error("Z_ERROR", multiple * 0),
        "split" : Error("Z_ERROR", multiple * 0),
        "shift" : Error("Z_ERROR", multiple * 0), 
        "shift_prop_to" : None,

        # Opetional when we were considering just pausing the syndrome extraction.
        "pause" : Error("Z_ERROR", multiple * 0),
        
        "four_ion_shift" : Error("Z_ERROR", multiple * p_idle_dephasing(t_4s, T2), multiple * m * t_4s, round_sig_fig( multiple * m * t_4s / (1 - multiple * m * t_4s), 5)),  
    }

    return helios_idle_during


''' zero_errors
    Set all errors to zero'''
def zero_errors():
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


        # Added repumping error from Honeywell in case it gets called:
        "error_per_repump" : Error("DEPOLARIZE1", 0) 
        
    }
    
    return errors


''' zero_idling
    Defines default idling gates and error rates, i.e. the operation and probability of that operation applied to qubits that are idling in a timestep that other qubits are experiencing the key operation'''
def zero_idling():

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
    
    p_relax = round_sig_fig(( p / (1 - p)), 6)
    p_leak = p

    errors = {

        "RZ" : Error("X_ERROR", p, p_leak, p_relax), 
        "RX" : Error("Z_ERROR", p, p_leak, p_relax),
        "H" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "CNOT" : Error("DEPOLARIZE2", p, p_leak, p_relax),
        "CZ" : Error("DEPOLARIZE2", p, p_leak, p_relax),
        "MZ" : Error("X_ERROR", p, p_leak, p_relax),
        "MX" : Error("Z_ERROR", p, p_leak, p_relax),


        # Qubit module errors -- None
        "shift" : Error("DEPOLARIZE1", 0),
        "shift_prop_to" : None, 

        # Additional for our architecture - None
        "shuttle" : Error("DEPOLARIZE1", 0), 
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),

        # Added repumping error from Honeywell in case it gets called:
        "error_per_repump" : Error("DEPOLARIZE1", p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 

    }
    
    return errors


def uniform_idling(p):


    p_relax = round_sig_fig(( p / (1 - p)), 6)
    p_leak = p

    idle_during = {
        "RZ" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "RX" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "H" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "CNOT" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "CZ" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "MZ" : Error("DEPOLARIZE1", p, p_leak, p_relax),
        "MX" : Error("DEPOLARIZE1", p, p_leak, p_relax),

        "shuttle" : Error("DEPOLARIZE1", 0),
        "merge" : Error("DEPOLARIZE1", 0),
        "split" : Error("DEPOLARIZE1", 0),
        "shift" : Error("DEPOLARIZE1", 0),

        "shift_prop_to" : 0,

        "pause" : Error("DEPOLARIZE1", 0),

    }

    return idle_during



''' round_sig_fig
Rounds a number to the given amount of significant figures'''
def round_sig_fig(x, sig):
    if x == 0:
        return 0
    else:
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)



''' p_idle_dephasing
If a qubit is idling for time t, this function returns what the probability of a Z-error occurring on it will be. It takes as inputs t, the time the qubit is idling for, and T2, the characteristic time for dephasing of an idling qubit. For example p_idle(t = 100e-6, T2 = 50s) = 1e-6 , indicating that if a qubit is idling for 100μs it will experience a Z-error with probability 1e-6. This function assumes idling is a dephasing noise channel with p = 0.5(1 - e^(-t/T_2))'''
def p_idle_dephasing(t, T2):
    p = 0.5 * (1 - np.exp(-t / T2))
    p_rounded = round_sig_fig(p, 3) # round to 3 sig fig
    return p_rounded

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


        # Added repumping error from Honeywell in case it gets called:
        "error_per_repump" : Error("DEPOLARIZE1", 2e-2 * p)  # For every repumping cycle there is a memory error on non-leaked qubits. This is not specified in the Quantinuum Helios paper but is instead based on this earlier paper from them (then Honeywell):  10.1103/PhysRevLett.124.170501 
    }
    
    return errors

