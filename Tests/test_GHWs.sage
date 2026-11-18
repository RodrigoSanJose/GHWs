# Required imports

import time

# Load the package GHWs

from GHWs import *

# Formula for the GHWs of Reed-Muller codes

def GHWRM(q, d, m, r):
    dt = m * (q - 1) - d - 1
    A = IntegerListsLex(min_sum=dt + 1, max_part=q - 1, length=m)
    A = A.list()
    A.reverse()
    try:
        j = A[r - 1]
    except:
        raise Exception("r is greater than the dimension")
    temp=1
    for y in srange(m):
        temp= temp + j[y]*q^(m - y -1)
    return temp

start = time.time()

# We test the functions hierarchy and GHW with the GHWs
# of some binary Reed-Muller codes

q = 2
mlimit = 6
for m in srange(2, mlimit+1):
    for d in srange(floor((m - 1) / 2)): 
        C = codes.BinaryReedMullerCode(d, m)
        k = C.dimension()
        n = C.length()
        hie = [GHWRM(q, d, m, r) for r in srange(1, k + 1)]
        if not hie == [GHW(C, r) for r in srange(1, k + 1)] or not hie == hierarchy(C):
            raise Exception("Error with GHW or hierarchy in RM test")


# We test the functions rhierarchy and RGHW with the
# RGHWs of some q-ary Reed-Muller codes

C1 = LinearCode(Matrix(GF(5), [(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
 (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
 (0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4),
 (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1, 1, 1),
 (0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 0, 2, 4, 1, 3, 0, 3, 1, 4, 2, 0, 4, 3, 2, 1),
 (0, 1, 4, 4, 1, 0, 1, 4, 4, 1, 0, 1, 4, 4, 1, 0, 1, 4, 4, 1, 0, 1, 4, 4, 1)]))

C2 = LinearCode(Matrix(GF(5), [(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
 (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
 (0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4)]))

rhie = rhierarchy(C1, C2)
if not rhie == [15,19,22] or not rhie == [RGHW(C1, C2, r) for r in srange(1, C1.dimension() - C2.dimension() + 1)]:
    raise Exception("Error with the RGHW or rhierarchy in RM test")

# We test GHW, hierarchy and wei_duality with the
# GHWs of Reed-Solomon codes
QQ = [i for i in srange(2, 14) if i.is_prime_power()]
for q in QQ:
    K = GF(q)
    for d in srange(0, floor((q - 1) / 2) + 1): 
        C = codes.GeneralizedReedSolomonCode(K.list(), d + 1)
        k = d + 1
        n = q
        hie = [n - k + r for r in srange(1, k + 1)]
        if not hie == [GHW(C, r) for r in srange(1, k + 1)] or not hie == hierarchy(C):
            raise Exception("Error with GHW or hierarchy in RS test")
        hie_dual = [i for i in srange(1, n + 1) if n + 1 - i not in hie]
        if not hie_dual == wei_duality(hie):
            raise Exception("Error with wei_duality in RS test")
        
# We test GHW, and the functions related to cyclic codes,
# with the GHWs of some binary cyclic codes

C = codes.BCHCode(GF(2), 35, 7)
if not [GHW(C, r) for r in range(1, 4)] == [7, 14, 21]:
    raise Exception("Error with GHW in the cyclic test")
C = codes.BCHCode(GF(2), 63, 27)
if not [GHW(C, r) for r in range(1, 4)] == [27, 41 ,48]:
    raise Exception("Error with GHW in the cyclic test")

# We test higher_spectrum and rhigher_spectrum with
# a Reed-Solomon code

q = 7
K = GF(q)
n = q
k = 4
k2 = 2 
C = codes.GeneralizedReedSolomonCode(K.list(), k)
C2 = codes.GeneralizedReedSolomonCode(K.list(), k2)
spec, esp = higher_spectrum(C, subspace_list=True)
rspec, resp = rhigher_spectrum(C, C2, subspace_list=True)
if not spec == {0: {0: 1}, 1: {4: 35, 5: 63, 6: 168, 7: 134}, 2: {5: 21, 6: 357, 7: 2472}, 3: {6: 7, 7: 393}, 4: {7: 1}}:
    raise Exception("Error with higher_spectrum in the RS test")
if not rspec == {0: {0: 1}, 1: {4: 35, 5: 63, 6: 161, 7: 133}, 2: {5: 21, 6: 301, 7: 2079}}:
    raise Exception("Error with rhigher_spectrum in the RS test")

# Check the order of matrices when two subspaces have the same support size.
K = GF(2)
A = matrix(K, [[1, 0, 1, 0]])
B = matrix(K, [[0, 1, 0, 1]])
C = LinearCode(A.stack(B))
C2 = LinearCode(A + B)
expected = [{1: {2: 2, 4: 1}}, {1: {2: [A, B], 4: [A + B]}}]
for spectrum in (higher_spectrum, higher_spectrum_low_mem):
    if spectrum(C, R=[1], subspace_list=True) != expected:
        raise Exception(f"Error with the subspace list in {spectrum.__name__}")
expected_relative = [{1: {2: 2}}, {1: {2: [A, B]}}]
for spectrum in (rhigher_spectrum, rhigher_spectrum_low_mem):
    if spectrum(C, C2, R=[1], subspace_list=True) != expected_relative:
        raise Exception(f"Error with the subspace list in {spectrum.__name__}")

# We test the strategies of hierarchy (direct, or with the dual code and Wei
# duality) with 2k > n, 2k = n, a cyclic code, a degenerate code, a code with
# a degenerate dual, and a code that is the whole space on its support

from contextlib import redirect_stdout
from io import StringIO

K = GF(2)
H = codes.BCHCode(K, 7, 3).generator_matrix() # [7, 4, 3] cyclic Hamming code
R = codes.BinaryReedMullerCode(1, 3).generator_matrix() # [8, 4, 4] code

def auto_route(C, L=None): # Hierarchy given by 'auto', and whether it used the dual code
    out = StringIO()
    with redirect_stdout(out):
        hie = hierarchy(C, L=L, verbose=True)
    return hie, 'Wei duality' in out.getvalue()

# Each test: generator matrix, weight hierarchy, and whether 'auto' uses the dual code
tests = [(H, [3, 5, 6, 7], True), (H.augment(matrix(K, 4, 2)), [3, 5, 6, 7], True), (R, [4, 6, 7, 8], False),
         (R.block_sum(matrix(K, [[1]])), [1, 5, 7, 8, 9], True),
         (identity_matrix(K, 3).augment(matrix(K, 3, 1)), [1, 2, 3], True)]
for G, hie, dual in tests:
    C = LinearCode(G)
    if not auto_route(C) == (hie, dual) or not all(hierarchy(C, strategy=s) == hie for s in ['direct', 'dual']) or not hie == [GHW(C, r) for r in srange(1, C.dimension() + 1)]:
        raise Exception("Error with the strategies of hierarchy")
# With L given (here from a redundant generating set), 'auto' uses the direct route, although 2k > n for H
if not auto_route(LinearCode(H), information(H.stack(H))) == ([3, 5, 6, 7], False):
    raise Exception("Error with the strategies of hierarchy")
for args in [{'L': information(H), 'strategy': 'dual'}, {'strategy': 'other'}, {'L': [[[0, 1, 2]], [H[:3]], [0]]}]:
    try:
        hierarchy(LinearCode(H), **args)
    except Exception:
        continue
    raise Exception("Error with the arguments of hierarchy")

end = time.time()
print('Total time:', end - start)