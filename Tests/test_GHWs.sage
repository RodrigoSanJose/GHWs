# Required imports

import time
from pathlib import Path

# Load the package GHWs

path1 = Path(__file__).resolve()
path2 = path1.parent.parent / 'GHWs' / 'GHWs.py'
load(path2)

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

end = time.time()
print('Total time:', end - start)