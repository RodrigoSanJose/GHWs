# Required imports

from pathlib import Path

# Load the package GHWs

path1 = Path(__file__).resolve()
path2 = path1.parent.parent / 'GHWs' / 'GHWs.py'
load(path2)

# Naive algorithm

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

r = 2
Q = [2, 3, 4, 5, 7, 8]
C = codes.HammingCode(GF(2), 5).dual_code()
G = C.generator_matrix()
for q in Q:
    C = LinearCode(matrix(GF(q), G))
    print('q =', q, 'GHWs:', timeit('GHW_low_mem(C,r)', precision=10, number=25, repeat=5))
    print('q =', q, 'Naive method:', timeit('naive_GHW_low_mem(C,r)', precision=10, number=25, repeat=5))
    print('-'*60)