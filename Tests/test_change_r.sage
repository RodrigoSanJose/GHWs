# Required imports

import time
from itertools import combinations
from sage.coding.reed_muller_code import QAryReedMullerCode
from pathlib import Path

# Load the package GHWs

path1 = Path(__file__).resolve()
path2 = path1.parent.parent / 'GHWs' / 'GHWs.py'
load(path2)

# Formula for the GHWs of Reed-Muller codes

def GHWRM(q, d, m, r):
    dt = m*(q-1) - d - 1
    A = IntegerListsLex(min_sum=dt + 1,max_part=q - 1,length=m)
    A = A.list()
    A.reverse()
    try:
        j = A[r - 1]
    except:
        raise Exception("r is greater than the dimension")
    temp = 1
    for y in srange(m):
        temp = temp + j[y]*q^(m - y - 1)
    return temp

# Naive algorithm

def naive_GHW(C, r): # Standard version
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception("Invalid value of r")
    w = r
    ghwub = n - k + r
    while w <= k: 
        rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
        it = 0
        for y in combinations(range(k), w): # All possible supports of weight w
            for mat in rrefs:
                it = it +1
                MM = []
                for t in range(k):
                    if t in y:
                        MM.append(mat.column(y.index(t))) 
                    else:
                        MM.append([0 for z in range(r)])
                ghwub = min(ghwub, len(matrix_supp(matrix(MM).transpose()*G)))
        w = w + 1
    return ghwub

def naive_GHW_low_mem(C, r): # Low memory version
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception("Invalid value of r")
    w = r
    ghwub = n - k + r
    while w <= k: 
        y = range(w) # We start with support {1,...,w}
        for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
            # position, and we choose r-1 more pivots
            comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
            ss = [y[0]]+list(s) # Pivots
            weights_columns = []
            for i in comp:
                j = 1
                while i - j not in ss:
                    j = j + 1
                ind = ss.index(i - j) 
                weights_columns.append(ind + 1) # The columns can have
                # different weights depending on their position relative
                # to the pivots
            # The list of all possible non-pivot columns:
            fqcols = cartesian_product([colwt(j, r, K) for j in weights_columns])
            for i in fqcols:
                for sup in combinations(range(k), w): # All possible supports of weight w
                    MM = []
                    for t in range(k):
                        if t in sup:
                            if sup.index(t) in ss: # Pivots
                                MM.append(standard(ss.index(sup.index(t)), r, K))
                            else: # Non pivots
                                MM.append(i[comp.index(sup.index(t))])
                        else:
                            MM.append([0 for z in range(r)])
                    ghwub = min(ghwub, len(matrix_supp(matrix(MM).transpose() * G)))
        w = w + 1
    return ghwub

# Tests

Q = [5, 7]
R = [2, 3, 4, 5]
for q in Q:
    for r in R:
        C = QAryReedMullerCode(GF(q), 2, 2)
        start = time.time()
        g1 = GHW_low_mem(C, r)
        end = time.time()
        print('q =', q, 'r =', r, 'GHWs:', end - start)
        g2 = naive_GHW_low_mem(C, r) 
        end2 = time.time()
        if g1 != g2 or g1 != GHWRM(q,2,2,r):
            raise Exception('Error in the computation of the GHWs')
        print('q =', q, 'r =', r, 'Naive method:', end2 - end)
        print('-'*60)
        if q == 8 and r == 2:
            break 
