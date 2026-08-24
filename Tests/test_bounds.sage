# Load the package GHWs

from GHWs import *

# We test lower and upper bounds for GHWs with a Reed-Solomon code

K = GF(7)
C = codes.GeneralizedReedSolomonCode(K.list(), 4)
L = information(C.generator_matrix())

# The second GHW is 5. Thus, the lower bound 5 and the upper bound 5 are
# attained, whereas the lower bound 4 and the upper bound 6 are not.
if not GHW_bound(C, 2, 5, L=L) or GHW_bound(C, 2, 4, L=L):
    raise Exception('Error with lower bounds in GHW_bound test')
if not GHW_bound(C, 2, 5, bound_type='upper', L=L) or GHW_bound(C, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in GHW_bound test')
if not GHW_bound_low_mem(C, 2, 5, L=L) or GHW_bound_low_mem(C, 2, 4, L=L):
    raise Exception('Error with lower bounds in GHW_bound_low_mem test')
if not GHW_bound_low_mem(C, 2, 5, bound_type='upper', L=L) or GHW_bound_low_mem(C, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in GHW_bound_low_mem test')

# We test an upper bound that is contradicted by the BCH bound

C = codes.BCHCode(GF(2), 15, 3)
if GHW_bound(C, 1, 2, bound_type='upper') or GHW_bound_low_mem(C, 1, 2, bound_type='upper'):
    raise Exception('Error with the BCH bound in GHW_bound test')

# We test lower and upper bounds for RGHWs

G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
G2 = matrix(GF(2), [G[-1]])
C = LinearCode(G)
C2 = LinearCode(G2)
L = information(C.generator_matrix())

# The second RGHW is 5. Thus, the lower bound 5 and the upper bound 5 are
# attained, whereas the lower bound 4 and the upper bound 6 are not.
if not RGHW_bound(C, C2, 2, 5, L=L) or RGHW_bound(C, C2, 2, 4, L=L):
    raise Exception('Error with lower bounds in RGHW_bound test')
if not RGHW_bound(C, C2, 2, 5, bound_type='upper', L=L) or RGHW_bound(C, C2, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in RGHW_bound test')
if not RGHW_bound_low_mem(C, C2, 2, 5, L=L) or RGHW_bound_low_mem(C, C2, 2, 4, L=L):
    raise Exception('Error with lower bounds in RGHW_bound_low_mem test')
if not RGHW_bound_low_mem(C, C2, 2, 5, bound_type='upper', L=L) or RGHW_bound_low_mem(C, C2, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in RGHW_bound_low_mem test')

# These cases require finding an actual witness.
if not GHW_bound(C, 2, 4) or not GHW_bound_low_mem(C, 2, 4):
    raise Exception('Error finding a witness in GHW_bound test')

if GHW_bound(C, 2, 5, bound_type='upper') or GHW_bound_low_mem(
        C, 2, 5, bound_type='upper'):
    raise Exception('Error rejecting a non-tight upper GHW bound')

if not RGHW_bound(C, C2, 1, 3) or not RGHW_bound_low_mem(C, C2, 1, 3):
    raise Exception('Error finding a witness in RGHW_bound test')

if RGHW_bound(C, C2, 1, 4, bound_type='upper') or RGHW_bound_low_mem(
        C, C2, 1, 4, bound_type='upper'):
    raise Exception('Error rejecting a non-tight upper RGHW bound')

if GHW(C, 2) != 4 or RGHW(C, C2, 1) != 3:
    raise Exception('Unexpected weights in bound test example')
