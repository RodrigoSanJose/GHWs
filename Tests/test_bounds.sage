# Required imports

import time

# Load the package GHWs

from GHWs import *

start = time.time()
# We test lower and upper bounds for GHWs with a Reed-Solomon code

K = GF(7)
C = codes.GeneralizedReedSolomonCode(K.list(), 4)
L = information(C.generator_matrix())

# The second GHW is 5. Hence, 4 and 5 are valid lower bounds, whereas 6 is
# not. Similarly, 5 and 6 are valid upper bounds, whereas 4 is not.
if not GHW_bound(C, 2, 4, L=L) or not GHW_bound(C, 2, 5, L=L) or GHW_bound(C, 2, 6, L=L):
    raise Exception('Error with lower bounds in GHW_bound test')
if GHW_bound(C, 2, 4, bound_type='upper', L=L) or not GHW_bound(C, 2, 5, bound_type='upper', L=L) or not GHW_bound(C, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in GHW_bound test')
if not GHW_bound_low_mem(C, 2, 4, L=L) or not GHW_bound_low_mem(C, 2, 5, L=L) or GHW_bound_low_mem(C, 2, 6, L=L):
    raise Exception('Error with lower bounds in GHW_bound_low_mem test')
if GHW_bound_low_mem(C, 2, 4, bound_type='upper', L=L) or not GHW_bound_low_mem(C, 2, 5, bound_type='upper', L=L) or not GHW_bound_low_mem(C, 2, 6, bound_type='upper', L=L):
    raise Exception('Error with upper bounds in GHW_bound_low_mem test')

# We test bounds that can be decided directly using the BCH bound

C = codes.BCHCode(GF(2), 15, 3)
if not GHW_bound(C, 1, 3) or not GHW_bound_low_mem(C, 1, 3):
    raise Exception('Error validating the BCH lower bound')
if GHW_bound(C, 1, 2, bound_type='upper') or GHW_bound_low_mem(C, 1, 2, bound_type='upper'):
    raise Exception('Error rejecting an upper bound contradicted by the BCH bound')

# We test cases in which an actual witness or counterexample has to be found

G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
G2 = matrix(GF(2), [G[-1]])
C = LinearCode(G)
C2 = LinearCode(G2)
L = information(C.generator_matrix())

# The first and second GHWs are 2 and 4, respectively. The invalid lower bound
# 5 and the valid upper bound 4 require finding a subcode with support 4. The
# invalid upper bound 3 is rejected once all the remaining subspaces have
# support at least 4. For r=1, the upper bound 3 is valid but not tight.
if not GHW_bound(C, 2, 3, L=L) or not GHW_bound(C, 2, 4, L=L) or GHW_bound(C, 2, 5, L=L):
    raise Exception('Error validating lower GHW bounds with witnesses')
if GHW_bound(C, 2, 3, bound_type='upper', L=L) or not GHW_bound(C, 2, 4, bound_type='upper', L=L):
    raise Exception('Error validating upper GHW bounds with witnesses')
if GHW_bound(C, 1, 3, L=L) or not GHW_bound(C, 1, 3, bound_type='upper', L=L):
    raise Exception('Error with non-tight GHW bounds')

if not GHW_bound_low_mem(C, 2, 3, L=L) or not GHW_bound_low_mem(C, 2, 4, L=L) or GHW_bound_low_mem(C, 2, 5, L=L):
    raise Exception('Error validating lower GHW bounds with witnesses in low memory')
if GHW_bound_low_mem(C, 2, 3, bound_type='upper', L=L) or not GHW_bound_low_mem(C, 2, 4, bound_type='upper', L=L):
    raise Exception('Error validating upper GHW bounds with witnesses in low memory')
if GHW_bound_low_mem(C, 1, 3, L=L) or not GHW_bound_low_mem(C, 1, 3, bound_type='upper', L=L):
    raise Exception('Error with non-tight GHW bounds in low memory')

# The first two RGHWs are 3 and 5. We test tight and non-tight valid bounds,
# as well as invalid bounds on both sides.
if not RGHW_bound(C, C2, 1, 2, L=L) or not RGHW_bound(C, C2, 1, 3, L=L) or RGHW_bound(C, C2, 1, 4, L=L):
    raise Exception('Error validating lower RGHW bounds')
if RGHW_bound(C, C2, 1, 2, bound_type='upper', L=L) or not RGHW_bound(C, C2, 1, 3, bound_type='upper', L=L) or not RGHW_bound(C, C2, 1, 4, bound_type='upper', L=L):
    raise Exception('Error validating upper RGHW bounds')
if not RGHW_bound(C, C2, 2, 4, L=L) or not RGHW_bound(C, C2, 2, 5, L=L) or RGHW_bound(C, C2, 2, 4, bound_type='upper', L=L):
    raise Exception('Error validating bounds for the second RGHW')

if not RGHW_bound_low_mem(C, C2, 1, 2, L=L) or not RGHW_bound_low_mem(C, C2, 1, 3, L=L) or RGHW_bound_low_mem(C, C2, 1, 4, L=L):
    raise Exception('Error validating lower RGHW bounds in low memory')
if RGHW_bound_low_mem(C, C2, 1, 2, bound_type='upper', L=L) or not RGHW_bound_low_mem(C, C2, 1, 3, bound_type='upper', L=L) or not RGHW_bound_low_mem(C, C2, 1, 4, bound_type='upper', L=L):
    raise Exception('Error validating upper RGHW bounds in low memory')
if not RGHW_bound_low_mem(C, C2, 2, 4, L=L) or not RGHW_bound_low_mem(C, C2, 2, 5, L=L) or RGHW_bound_low_mem(C, C2, 2, 4, bound_type='upper', L=L):
    raise Exception('Error validating bounds for the second RGHW in low memory')

# Sanity checks for the weights used above
if GHW(C, 1, L=L) != 2 or GHW(C, 2, L=L) != 4:
    raise Exception('Unexpected GHWs in bound test example')
if RGHW(C, C2, 1, L=L) != 3 or RGHW(C, C2, 2, L=L) != 5:
    raise Exception('Unexpected RGHWs in bound test example')

end = time.time()
print('Total time:', end - start)
