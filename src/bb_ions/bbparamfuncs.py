''' bbparamfuncs
Once the parity check matrices of a BB code have been constructed using the functions in bbfuncs.py, it is useful to find the paramaters of the code.
I.e. how many logical qubits it contains, what each logical qubit's logical operators (X_L and Z_L) are and the distance of the code.

Note on code distance:
Code distance is the smallest number of single-qubit operations required to cause a logical error. In the case of this being different for different
logical qubits in the BB code, the smallest distance is given. Note however that the distance returned from the bposd simulations is a 'maximum' distance.
That is, while the smallest layout of single-qubit operations that it found that caused a logical error might be, for example, d = 5, it can not guarantee
that a layout with a smaller number, for example d = 4, doesn't exist. Hence d = 5 is the maximum possible distance that the code could have. The more 
trials that are run the more confident we become that this is its actual distance however it's not certain.
'''

import numpy as np
import bposd
from bposd import css
from bposd import css_decode_sim
import contextlib
from autqec.utils.qec import *
from .bbfuncs import *
import os
import sys
import random
from collections import defaultdict

# # bbq for tanner graph (use code.tanner_graph) but is not always the best layout (can have unnecessarily long connections) and not the same qubit coordinates as my stim circuits - commenting out for now but can be put back in (along with the commented out @property below) for a quick visualisation of BB codes, though circuit.diagram("timeslice-svg") is also good as basically reveals what the tanner graph is for my stim circuits)
# import bbq # install by changing to bbqudit directory and pip install e . (the pip install bbqudit version does not have Monomial etc.)
# from bbq.polynomial import Monomial
# from bbq.bbq_code import BivariateBicycle
# from bbq.field import Field
# import matplotlib.pyplot as plt

class Code:
    def __init__(self, l, m, Aij, Bij, ATij, BTij, Hx, Hz, Lx, Lz, d_max, n, k, Junion, JTunion, As_is, Bs_is, ATs_is, BTs_is, A, B):
        self.l = l
        self.m = m
        self.Aij = Aij
        self.Bij = Bij
        self.ATij = ATij
        self.BTij = BTij
        self.Hx = Hx
        self.Hz = Hz
        self.Lx = Lx
        self.Lz = Lz
        self.d_max = d_max
        self.n = n
        self.k = k
        self.Junion = Junion
        self.JTunion = JTunion
        self.As_is = As_is  # " A's i's -- i.e. a dictionary that you can put a jvalue into and it returns all the i values that go with that j value"
        self.Bs_is = Bs_is
        self.ATs_is = ATs_is
        self.BTs_is = BTs_is
        self.A = A
        self.B = B

    # @property
    # def tanner_graph(self):
    #     if hasattr(self, "_tanner_graph_done"):
    #         return

    #     field = Field(2)
    #     x = Monomial(field, 'x')
    #     y = Monomial(field, 'y')

    #     A = sum((x**i * y**j for (i, j) in self.Aij))
    #     B = sum((x**i * y**j for (i, j) in self.Bij))

    #     bb = BivariateBicycle(A, B, self.l, self.m, 1)
    #     bb.draw()
    #     plt.title(r"Tanner graph (dif. coords. to stim circuit)")

    #     self._tanner_graph_done = True


''' Functions defining BB codes'''


''' Some small codes I found'''
def bb6_6_4_2(): # Iceberg [[2m, 2m - 2, 2]] code. [[6, 4, 2]]. Its parity-check matrices are all 1's. Can be expressed as a BB code with the following parameters:
    l = 3
    m = 1
    Aij= [(0, 0), (1, 0), (2, 0)]
    Bij =[(0, 0), (1, 0), (2, 0)]
    d = 2
    code = get_code_params(l, m, Aij, Bij, d)
    return code

iceberg_code = bb6_6_4_2

def bb8_8_6_2(): #Iceberg 
    l = 2
    m = 2
    Aij = [[0, 0], [0, 1], [1, 0], [1, 1]]
    Bij = [[0, 0], [0, 1], [1, 0], [1, 1]]
    d = 2 
    code = get_code_params(l, m, Aij, Bij, d)
    return code 

def bb5_12_4_2():
    l = 3
    m = 2
    Aij = [[0, 0], [0, 1]]
    Bij = [[0, 0], [1, 1], [2, 0]]
    d = 2
    code = get_code_params(l, m, Aij, Bij, d)
    return code

def bb6_12_4_2():
    l= 3
    m = 2
    Aij= [(0, 0), (1, 0), (2, 0)]
    Bij =[(0, 0), (1, 1), (2, 0)]
    d = 2
    code = get_code_params(l, m, Aij, Bij, d)
    return code

def bb8_12_4_2():
    l = 3
    m = 2
    Aij = [[0, 0], [0, 1], [1, 0], [2, 1]]
    Bij = [[0, 0], [0, 1], [1, 0], [1, 1]]
    return get_code_params(l, m, Aij, Bij)

# def bb6_12_8_2():  ## Actually this is just three copies of the [[4, 2, 2]] code.
#     l= 3
#     m = 2
#     Aij= [(0, 0), (1, 0), (2, 0)]
#     Bij =[(0, 1), (1, 1), (2, 1)]
#     d = 2
#     code = get_code_params(l, m, Aij, Bij, d)
#     return code

# def bb8_12_8_2(): # Is just two copies of a BB8 [[6, 4, 2]] code
#     l = 3
#     m = 2
#     Aij = [[0, 0], [0, 1], [2, 0], [2, 1]]
#     Bij = [[0, 0], [0, 1], [1, 0], [1, 1]]
#     d = 2
#     code = get_code_params(l, m, Aij, Bij, d)
#     return code


''' Larger BB codes from papers'''

def bb5_18_4_3():
    # [[18, 4, 3]] BB5 code from IonQ experimental realisation of many qLDPC codes []
    l = 3
    m = 3
    Aij = [(0, 0), (1, 0)]
    Bij = [(0, 0), (0, 1), (1, 2)]
    d = 3
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb6_18_4_3():
    # [[18, 4, 3]] BB6 code from relayBP github: https://tinyurl.com/5a7urtbf
    l = 3 # only other options are l = 9 and m = 1 or vice versa which doesn't give logical qubits so must be  l = m = 3
    m = 3 
    # A = x + 1 + y
    # B = x + 1 + xy^2
    Aij = [(1,0), (0,0), (0,1)]
    Bij = [(1,0), (0,0), (1,2)]
    d=3 # smallest d I found was 3
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb4_18_2_3():
    # Small code example from Tripier ... Delfosse IonQ [2604.19481], Fig. 13
    l = 3 
    m = 3
    Aij = [(0,0), (1,0)]
    Bij = [(0,0), (1,1)]
    d = 3
    code = get_code_params(l,m,Aij,Bij,d)
    return code


def bb6_18_4_4():
    # [[18, 4, 4]] BB6 (weight-6 stabilisers) code from Wang et al. [2505.09684]
    l = 3
    m = 3
    # A = x + y^0 + y^2
    # B = y + x^0 + x^2
    Aij = [(1,0), (0,0), (0,2)]
    Bij = [(0,1), (0,0), (2,0)]
    d = 4
    code = get_code_params(l, m, Aij, Bij, d)

    return code

def bb8_18_6_3():
    l = 3 
    m = 3
    Aij = [(0, 0), (0, 1), (2, 1), (2, 2)]
    Bij = [(0, 2), (1, 0), (1, 1), (2, 2)]
    d = 3
    code = get_code_params(l, m, Aij, Bij, d)
    return code



'''Some found BB6 codes that might have been those compared to in the long chain paper [2503.22071]'''

'''bb6_30_4_4_code'''
def bb6_30_4_4():
    # Equal best of six [[30,4,4]] and [[30,4,6]] codes found using a search in "find_codes" directory
    l = 5
    m = 3
    Aij = [(0, 0), (0, 1), (1, 2)]
    Bij = [(0, 1), (0, 2), (2, 0)]
    d = 4
    code = get_code_params(l, m, Aij, Bij, 4)
    return code



'''bb6_48_4_6_code'''
def bb6_48_4_6():
    ## Best of three of the [[48,4,6]] codes found using a search in "find_codes" directory 
    l = 6
    m = 4 
    Aij = [(0, 0), (1, 0), (2, 1)]
    Bij = [(0, 1), (1, 0), (2, 1)]
    d = 6
    code = get_code_params(l, m, Aij, Bij, d)

    return code


''' BB5 codes from BB5 long chain paper [2503.22071] '''

def bb5_30_4_5():
    
    # [[30, 4, 5]] BB5 (weight-5 stabilisers) code from Ye Delfosse long chain [2503.22071] Table II
    l = 5
    m = 3
    # A = x^0 + x
    # B = x^0 + y + x^2*y^2
    Aij = [(0, 0), (1, 0)]          # the powers (i, j) of each term x^i * y^j in A
    Bij = [(0, 0), (0, 1), (2, 2)]
    d = 5
    code = get_code_params(l, m, Aij, Bij, d)
    
    return code


def bb5_48_4_7():
    # [[48, 4, 7]] BB5 code from Ye Delfosse long chain [2503.22071] Table II
    l = 8
    m = 3
    # A = x^0 + x
    # B = x^0 + y + x^3 * y^2
    Aij = [(0, 0), (1, 0)]
    Bij = [(0, 0), (0, 1), (3, 2)]
    d = 7

    code = get_code_params(l, m, Aij, Bij, d)

    return code



''' bb5 codes perform better'''
def bb5_60_4_8():
    l = 15
    m = 2
    # Aij = [[2, 0], [8, 1]]
    # Bij = [[0, 1], [4, 1], [14, 1]]
    Aij = [[2, 1], [8, 0]] # minusing 1 off of every y term just permutes the rows. Makes for nice polynomials.
    Bij = [[0, 0], [4, 0], [14, 0]]  
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code

def interchanged_bb5_60_4_8():
    l = 2
    m = 15
    Aij = [[1, 2], [0, 8]] # minusing 1 off of every y term just permutes the rows. Makes for nice polynomials.
    Bij = [[0, 0], [0, 4], [0, 14]]  
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code




''' Some found n≈60 BB6 codes that could also be run on Helios because they take n data qubits and n/2 check qubits so if (n + n/2) ≤ 98 Helios qubits'''
def bb8_48_6_8():
    l = 8
    m = 3
    Aij = [(2, 1), (3, 2), (5, 2), (7, 1)]
    Bij = [(1, 2), (2, 0), (5, 0), (5, 2)]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb6_56_6_8():
    l = 7
    m = 4
    Aij = [[0, 0], [2, 0], [3, 1]]
    Bij = [[3, 0], [5, 3], [6, 2]]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb6_60_8_6():
    l = 6
    m = 5  # NOTE this will be simulated as per our proposed architecture: one operation zone per module (so 5 in this case) but if pretending like this is being run on actual Helios, Helios only has 4 operation zones with 2q gates enabled so this is allowing for one more than that (won't affect pL very much as compared to doing completely sequential operations (imagining all disjoint operations can be done in a single time step) doesn't affect pL very much).
    Aij = [[0, 0], [0, 2], [1, 1]]
    Bij = [[0, 1], [2, 4], [3, 2]]
    d = 6
    code = get_code_params(l, m, Aij, Bij, d)
    return code


'''bb8 is not so great'''
def bb8_64_12_8():
    l = 8
    m = 4
    Aij = [(0, 2), (0, 3), (1, 1), (3, 2)]
    Bij = [(3, 3), (5, 2), (6, 1), (6, 2)]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code



''' End of some found codes around n=60'''



'''bb6_72_12_6'''
def bb6_72_12_6():
    # [[72, 12, 6]] BB6 code from OG BB paper [2308.07915] Table II
    l = 6
    m = 6
    # A = x^3 + y + y^2
    # B = y^3 + x + x^2
    Aij = [(3, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (1, 0), (2, 0)]
    d = 6
    code = get_code_params(l, m, Aij, Bij, d)
    return code

'''bb6_90_8_10_code'''
def bb6_90_8_10():
    # [[90, 8, 10]] BB6 code from OG BB paper [2308.07915] Table III
    l = 15
    m = 3
    # A = x^9 + y + y^2
    # B = x^0 + x^2 + x^7
    Aij = [(9, 0), (0, 1), (0, 2)]
    Bij = [(0, 0), (2, 0), (7, 0)]
    d = 10
    code = get_code_params(l, m, Aij, Bij, d)
    return code

def interchanged_bb6_90_8_10():
    # [[90, 8, 10]] BB6 code from OG BB paper [2308.07915] Table III
    l = 3
    m = 15
    Aij = [(0, 9), (1, 0), (2, 0)]
    Bij = [(0, 0), (0, 2), (0, 7)]
    d = 10
    code = get_code_params(l, m, Aij, Bij, d)
    return code




'''bb6_108_code'''
def bb6_108():
    # [[108, 8, 10]] BB6 code from OG BB paper [2308.07915] Table III
    l = 9
    m = 6
    # A = x^3 + y + y^2
    # B = y^3 + x + x^2
    Aij = [(3, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (1, 0), (2, 0)]
    d = 10
    code = get_code_params(l, m, Aij, Bij, d)
    return code

'''bb6_108_code'''
def interchanged_bb6_108():
    # [[108, 8, 10]] BB6 code from OG BB paper [2308.07915] Table III
    l = 6
    m = 9
    Aij = [(0, 3), (1, 0), (2, 0)]
    Bij = [(3, 0), (0, 1), (0, 2)]
    d = 10
    code = get_code_params(l, m, Aij, Bij, d)
    return code


'''bb5_120_8_8_code'''
def bb5_120_8_8():
    ## [[120, 8, 8]] BB5 found from search in "find_codes" directory
    l = 6
    m = 10
    Aij = [(0, 0), (0, 1)]
    Bij = [(0, 0), (2, 0), (4, 4)]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code


'''gross_code'''
def gross_code():
    # [[144, 12, 12]] BB6 'gross' code from OG BB paper [2308.07915] Table III
    l = 12
    m = 6
    # A = x^3 + y + y^2
    # B = y^3 + x + x^2
    Aij = [(3, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (1, 0), (2, 0)]
    d = 12

    code = get_code_params(l, m, Aij, Bij, d)

    return code

'''gross_code'''
def interchanged_gross_code():
    # [[144, 12, 12]] BB6 'gross' code from OG BB paper [2308.07915] Table III
    l = 6
    m = 12
    # A = x^3 + y + y^2
    # B = y^3 + x + x^2
    Aij = [(0, 3), (1, 0), (2, 0)]
    Bij = [(3, 0), (0, 1), (0, 2)]
    d = 12

    code = get_code_params(l, m, Aij, Bij, d)

    return code

'''bb6_210_code'''
def bb6_210_code():
    # A found 2-undercover of the found [[420, 16, 18]] code
    entry = {'nkd': [210, 8, 18], 'l': 7, 'm': 15, 'Aij': [(0, 0), (1, 2), (1, 8)], 'Bij': [(2, 13), (3, 11), (5, 4)]}
    return get_code_params(entry['l'], entry['m'], entry['Aij'], entry['Bij'], entry['nkd'][2])
bb6_210_8_18 = bb6_210_code

'''two_gross_code'''
def two_gross_code():
    # [[288, 12, 18]] BB6 'two gross' code from OG BB paper [2308.07915] Table III
    l = 12
    m = 12
    # A = x^3 + y^2 + y^7
    # B = y^3 + x^1 + x^2
    Aij = [(3, 0), (0, 2), (0, 7)]
    Bij = [(0, 3), (1, 0), (2, 0)]
    d = 18

    code = get_code_params(l, m, Aij, Bij, d)

    return code



def bb6_288_16_12():
    # 2-cover of gross code (can run ''find_codes/scripts/find_h_cover_codes.py gross_code 2 16 12''). (Identical to gross code except m = 2m_gross). Was found first using an LLM in [2606.02418v1] but really they should have searched h-covers first as it's much faster and less costly. I have simulated this code and sadly its performance isn't very good.
    l = 12
    m = 12
    Aij = [(3, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (1, 0), (2, 0)]
    d = 12
    code = get_code_params(l, m, Aij, Bij, d)

    # Is a multiple of the gross code? This is just a 2-cover of the gross code

    # l = 12 , m = 6,
    # A = x^3 + y + y^2
    # B = y^3 + x + x^2
    return code

two_cover_gross = bb6_288_16_12



def sum_gross_codes():
    # [[288, 24, 12]]
    # Is a direct sum of two gross codes -- from evolutionary discovery of BB codes paper [2606.02418v1]
    l = 12
    m = 12
    Aij = [(6, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (2, 0), (4, 0)]
    d = 12
    code = get_code_paramas(l, m, Aij, Bij, d)
    return code


'''bb6_360_code'''
def bb6_360_code():
    # [[360, 12, ≤24]] BB6 code from OG BB paper [2308.07915] Table III
    l = 30
    m = 6
    # A = x^9 + y^1 + y^2
    # B = y^3 + x^25 + x^26
    Aij = [(9, 0), (0, 1), (0, 2)]
    Bij = [(0, 3), (25, 0), (26, 0)]
    d = 24

    code = get_code_params(l, m, Aij, Bij, d)

    return code

def interchanged_bb6_360_code():
    # [[36x0, 12, ≤24]] BB6 code from OG BB paper [2308.07915] Table III BUT WITH l AND m (AND POWERS OF X AND Y) INTERCHANGED 
    # This makes a difference in performance as we are simulating m data-qubit modules of size 2*l and m operation zones when doing helios noise.
    l = 6
    m = 30
    Aij = [(0, 9), (1, 0), (2, 0)]
    Bij = [(3, 0), (0, 25), (0, 26)]
    d = 24
    code = get_code_params(l, m, Aij, Bij, d)
    return code


''' Instead of h-cover, searching for an undercover of the 756 code found:'''
def bb6_378_16_18():
    l= 21
    m = 9
    Aij = [(3, 0), (0, 10), (0, 17)]
    Bij = [(0, 5), (3, 0), (19, 0)]
    d_max = 18
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code
undercover_code = bb6_378_16_18


def interchanged_bb6_378_16_18():
    l= 9
    m = 21
    Aij = [(0, 3), (10, 0), (17, 0)]
    Bij = [(5, 0), (0, 3), (0, 19)]
    d_max = 18
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code
interchanged_undercover_code = interchanged_bb6_378_16_18




'''Codes from Abhishek Rajput and Ben Symons' paper 2511.13560'''

def bb8_64_14_8():
    # From 2511.13560
    l = 8
    m = 4
    # A = x y^3 + 1 + x^6 + x^3 y^2
    Aij = [(1,3),(0,0),(6,0),(3,2)]
    # B = x^6 y + x^4 y + x^3 + x^5 y
    Bij = [(6,1),(4,1),(3,0),(5,1)]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb8_72_14_8():
    # From 2511.13560
    l = 6
    m = 6
    # A = y^4 + x^5 y^4 + x^3 + x^5 y^3
    Aij = [(0,4),(5,4),(3,0),(5,3)]
    # B = y^5 + x^2 y + x^5 y^5 + x^3 y^4
    Bij = [(0,5),(2,1),(5,5),(3,4)]
    d = 8
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb8_128_14_12():
    # From 2511.13560
    l = 8
    m = 8
    # A = x y^3 + y^4 + x^6 y^4 + x^3 y^6
    Aij = [(1,3),(0,4),(6,4),(3,6)]
    # B = x^6 y + x^4 y^5 + x^3 + x^5 y^5
    Bij = [(6,1),(4,5),(3,0),(5,5)]
    d = 12
    code = get_code_params(l, m, Aij, Bij, d)
    return code

def fortnight_code():
    # From 2511.13560 - BB8 code
    l = 12
    m = 6
    # A = x^6 y^4 + x^5 y^4 + x^3 + x^11 y^3
    Aij = [(6,4),(5,4),(3,0),(11,3)]
    # B = y^5 + x^8 y + x^5 y^5 + x^9 y^4
    Bij = [(0,5),(8,1),(5,5),(9,4)]
    d = 14
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb8_192_14_16():
    # From 2511.13560 -- called 192_14_20 code but came back as d_max = 16 
    l = 8
    m = 12
    # A = x y^3 + y^8 + x^6 y^4 + x^3 y^2
    Aij = [(1,3),(0,8),(6,4),(3,2)]
    # B = x^6 y + x^4 y^5 + x^3 + x^5 y
    Bij = [(6,1),(4,5),(3,0),(5,1)]
    d = 16
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb8_216_14_20():
    # From 2511.13560 -- actually returns code distance of d_max ≤ 16.
    l = 18
    m = 6
    # A = x^6 y^4 + x^5 y^4 + x^15 + x^11 y^3
    Aij = [(6,4),(5,4),(15,0),(11,3)]
    # B = y^5 + x^8 y + x^11 y^5 + x^9 y^4
    Bij = [(0,5),(8,1),(11,5),(9,4)]
    d = 20
    code = get_code_params(l, m, Aij, Bij, d)
    return code


def bb6_248_10_14():
    # From 2511.13560 -- had d = 18 in paper but actually returns code distance of d_max ≤ 14.
    l = 31
    m = 4
    Aij = [(0,0),  (6,0), (27,0)]
    Bij = [(0,2), (15,3), (24,0)]
    d = 14
    code = get_code_params(l, m, Aij, Bij)
    return code

def bb6_648_code():
    l = 18
    m = 18
    d_max = 36
    Aij = [(3,0), (0, 13), (12, 2)]
    Bij = [(0,3), (7, 12), (14, 6)]
    return get_code_params(l, m, Aij, Bij, d_max)

'''Bigger BB8 codes subsequently found by Ben Symons'''

def bb8_288_14_20():
    l = 12
    m = 12
    Aij = [(6,10), (5,4), (3,6), (11,9)]
    Bij = [(0,5), (2,1), (11,5), (9,4)]
    d_max = 20
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code

def bb8_360_14_30():
    l = 30
    m = 6
    d_max = 30
    Aij = [(6,4), (29,4), (9,0), (5,3)]
    Bij = [(6,5), (20,1), (23,5), (21,4)]
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code


def four_twenty_code():
    # Found in this work.
    entries = [   # Statistically had pL all overlapping (four_twenty.csv and four_twenty_jupiter.csv) error bars but have left the 'best' one uncommented
        # {"nkd": [420, 16, 18], "l": 14, "m": 15, "Aij": [[0, 0], [4, 8], [6, 6]],  "Bij": [[5, 0], [5, 9], [7, 2]]}, 
        # {"nkd": [420, 16, 18], "l": 14, "m": 15, "Aij": [[0, 0], [0, 3], [8, 14]], "Bij": [[1, 10], [5, 13], [7, 9]]}, 
        # {"nkd": [420, 16, 18], "l": 14, "m": 15, "Aij": [[0, 0], [2, 2], [2, 9]],  "Bij": [[5, 3], [9, 0], [11, 4]]},
        {"nkd": [420, 16, 18], "l": 14, "m": 15, "Aij": [[0, 0], [8, 2], [8, 8]],  "Bij": [[2, 13], [10, 11], [12, 4]]},
        # {"nkd": [420, 16, 18], "l": 10, "m": 21, "Aij": [[0, 0], [2, 17], [8, 11]], "Bij": [[0, 15], [2, 15], [6, 11]]}
        ]
    entry = entries[0]
    return get_code_params(entry['l'], entry['m'], entry['Aij'], entry['Bij'], entry['nkd'][2])
bb6_420_code = four_twenty_code


'''bb6_756_code'''
def bb6_756_code():
    # [[756, 16, ≤34]] BB6 code from OG BB paper [2308.07915] Table III
    # Takes about a minute to find all the logical operators of this code
    l = 21
    m = 18
    Aij = [(3,0),(0,10),(0,17)]
    Bij = [(0,5),(3,0),(19,0)]
    d_max = 34
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code

def interchanged_756_code():
    # [[756, 16, ≤34]] BB6 code from OG BB paper [2308.07915] Table III
    # Takes about a minute to find all the logical operators of this code
    l = 18
    m = 21
    Aij = [(0,3),(10,0),(17,0)]
    Bij = [(5,0),(0,3),(0,19)]
    d_max = 34
    code = get_code_params(l, m, Aij, Bij, d_max)
    return code



''' suppress_stdout
For supressing the output from css_decode_sim'''
@contextlib.contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout




''' find_d_max
Given the X and Z parity-check matrices of a css code, this function finds the maximum distance it could have. I.e. it runs simulations with bposd decoder to find the distance of the code. We call this maximum distance as there could be a logical operator with a smaller weight that would cause a logical error (i.e., a smaller distance) but the simulations didn't see it. Increase target_runs for more certainty in d_max'''
def find_d_max(Hx, Hz, target_runs = 500):
    
    osd_options={
    'target_runs': target_runs, 'xyz_error_bias': [1, 1, 1], 'bp_method': "minimum_sum", 'ms_scaling_factor': 0.05, 'osd_method': "osd_cs", 'osd_order': 4, 'channel_update': None, 'seed': int(round(100*random.random())), 'max_iter': 9, 'error_bar_precision_cutoff': 1e-6, 'tqdm_disable' : 1
    }

    error_rate = 0.1

    with suppress_stdout():
        the_code = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = error_rate, **osd_options) 

    d_max = the_code.min_logical_weight 

    return d_max



''' get_code_params
Using the l, m, A indices (Aij) and B incdices (Bij) of a bicycle bivariate code this function constructs a code class which contains its parity check matrices, logical operators, d_max (could be a smaller distance but the simulations found d_max) n, k l, m, A and B.
Note that Aij and Bij are inserted are tuples (i, j) for each of the terms x^i y^j in A and B. 
E.g. if
A = x^0 + x
B = x^0 + y + x^2*y^2   ([[30, 4, 5]] code from Table II of [2503.22071])
Then
Aij = [(0, 0), (1, 0)]
Bij = [(0, 0), (0, 1), (2, 2)]
The A and B matrices returned by this function are the actual matrices A, B'''
def get_code_params(l, m, Aij, Bij, d_max = None, target_runs = 1000):



    # Sorting indices into I(A), I(B), J(A), J(B):
    IA, JA = findIJ(Aij)
    IB, JB = findIJ(Bij)
    Junion = sorted(set(JA + JB))

    # Sorting indices into I(A^T), I(B^T), J(A^T), J(B^T):
    ATij = [(-i, -j) for i, j in Aij]
    BTij = [(-i, -j) for i, j in Bij]
    
    
    IAT, JAT = findIJ(ATij)
    IBT, JBT = findIJ(BTij)
    JTunion = sorted(set(JAT + JBT))


    # Num qubits:
    n = 2 * l * m

    Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)

    # Logical operators:
    Lx, Lz = autqec_logical_ops(Hx, Hz)

    # Num. logical qubits:
    k = len(Lx)

    if k == 0:
        print("This code has no logical qubits")
        d_max = None
        code = Code(l, m, Aij, Bij, ATij, BTij, Hx, Hz, Lx, Lz, d_max, n, k, Junion, JTunion)    
        
        return code

    if d_max == None:
        d_max = find_d_max(Hx, Hz, target_runs = target_runs) 

    As_is = group_by_j(Aij)
    Bs_is = group_by_j(Bij)
    ATs_is = group_by_j(ATij)
    BTs_is = group_by_j(BTij)

    A=' + '.join(f"x^{i}⋅y^{j}" for i, j in Aij)
    B=' + '.join(f"x^{i}⋅y^{j}" for i, j in Bij)


    code = Code(l, m, Aij, Bij, ATij, BTij, Hx, Hz, Lx, Lz, d_max, n, k, Junion, JTunion, As_is, Bs_is, ATs_is, BTs_is, A, B)

    return code




''' bposd_logical_ops --(Recommended to use autqec_logical_ops function instead to return a canonical set of logical ops)
Given the parity check matrices of a BB code, let's get the logical operators for each of its logical qubits.
This function does this using Joschka Roffe's bposd package.
Checking the anticommutation relations of the logical operators is done with anticommute_matrix = Lx @ Lz^T % 2. It checks each Lx against every Lz. The anticommute matrix shows that bposd returns a generating set of logical operators that is not canonical. I.e. the anticommute matrix has rank = k but is not I_k. Row-reducing the anticommute matrix until it is I_k (and performing the corresponding multiplication of logical operators) will give a canonical set where each logical operator anticommutes with its partner on the same logical qubit but commutes with all the others.
This is done within the autqec_logical_ops function so it may be used instead.
'''
def bposd_logical_operators(Hx, Hz):
    code = css.css_code(hx = Hx, hz = Hz)
    # n = code.N
    # k = code.K
    # d = code.D
    # print(f'[[{n}, {k}, {d}]]\n')

    # Look at logical ops:
    Lx = code.lx.toarray()
    Lz = code.lz.toarray()

    # Check anticommutations between logical operators (a 1 implies anticommutes)
    anticommute_matrix = Lx @ Lz.T % 2
    rank = np.linalg.matrix_rank(anticommute_matrix)
    # Note if each logical operator anticommutes with exactly one other then you get the identity matrix out
    # You might also get an equivalent binary matrix out, i.e. of the same rank. If its rank is equal to number of logical qubits then you can multiply logical operators together to eventually just get pairs that anticommute, i.e. turn the matrix into the identity matrix)
    assert(rank == code.K)
        # Correct anticommutation relations between logical operators

    return Lx, Lz

''' autqec_logical_ops
Given the parity check matrices of a BB code, let's get the logical operators for each of its logical qubits.
This function does this using Hasan Sayginel's autqec package.
Checking the anticommutation relations of the logical operators is done with anticommute_matrix = Lx @ Lz^T % 2.
It checks each Lx against every Lz.
We check that the anticommmute matrix is exactly the identity, implying we have a canonical set of logial operators, where each logical operator anticommutes with its partner on the same logical qubit but commutes with all the others.
Note: if n and k have already been found they can be fed to this function. Alternatively it will find n and k using Joschka Roffe's bposd package.
'''
def autqec_logical_ops(Hx, Hz, n = None, k = None):
    
    if n is None or k is None:
        code = css.css_code(hx = Hx, hz = Hz)
        n = code.N
        k = code.K
    
    ''' 
    First add zero matrices to convert parity check matrices to symplectic form for autqec. 
    e.g. XYZIZ in symplectic form 
    = [X part | Z part]
    = [X X I I I | I Z Z I Z]
    = [1 1 0 0 0 | 0 1 1 0 1]

    Now because our BB code is CSS the combined parity check matrix is simply
     [  Hx  |  0 
        0   | Hz  ]
    '''

    zeros = np.zeros_like(Hx)
    H_symp = np.array(np.vstack((np.hstack((Hx,zeros)),np.hstack((zeros,Hz)))),dtype=int)

    # Row reduce and find logical operators:
    H_symp_rref, _, transform_rows, transform_cols = rref_mod2(H_symp)
    H_symp_rref = H_symp_rref[~np.all(H_symp_rref == 0, axis=1)]
    H_symp_rref_og_basis = H_symp_rref@inv_mod2(transform_cols)
    assert H_symp_rref_og_basis.shape[0] == n - k
    assert H_symp_rref_og_basis.shape[1] == 2 * n

    G, LX_symplectic, LZ_symplectic, D = compute_standard_form(H_symp_rref_og_basis)

    # We have a CSS code so cut off the empty Z and X parts of the symplectic representations of Lx and Lz:
    Lx = LX_symplectic[:, :n]
    Lz = LZ_symplectic[:, n:]

    # Verify anticommutation of logical operators and that they're canonical:
    anticommute_matrix = (Lx @ Lz.T) % 2
    ident_matrix = np.eye(Lx.shape[0])
    assert(np.array_equal(anticommute_matrix, ident_matrix))

    return Lx, Lz



''' findIJ
Given a bivariate polynomial, for example 
P = x^0 + y + x^2*y^2  
And writing tuples (i, j) for each term x^i⋅y^j, for example
P = [(0, 0), (0, 1), (2, 1)]
This function returns
I(P) = list of powers of x in P
J(P) = list of powers of y in P'''
def findIJ(P):
    IP, JP = zip(*P)
    return IP, JP

''' group_by_j
Creates a dictionary containing all the i's for a particular j.
E.g. Aij = [(0, 0), (0, 1), (2, 1)]
As_is = group_by_j(Aij)    # "A's i's" 
As_is[1] = [0, 2]'''
def group_by_j(P):
    j_to_i = defaultdict(list)
    for i, j in P:
        j_to_i[j].append(i)
    return dict(j_to_i)