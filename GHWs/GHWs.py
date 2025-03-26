#################################################################################
# V1.0 20/03/25
# Author:
# Rodrigo San-José. Contact: rodrigo.san-jose@uva.es
# GitHub repository: https://https://github.com/RodrigoSanJose/GHW

# This package provides functions to compute the generalized Hamming weights (GHWs)
# of a linear code, as well as its relative generalized Hamming weights, using a
# Brouwer-Zimmermann type of algorithm. 

# ****************************************************************************
#       Copyright (C) 2025 Rodrigo San-José <rodrigo.san-jose@uva.es>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ****************************************************************************

# Required imports
from itertools import combinations

#################################################################################
# AUXILIARY FUNCTIONS
#################################################################################

def vecwt(w, n, K): 
    r"""
    Computes all the vectors of weight w and length n over the finite field K

    OUTPUT: 
    
    A list of such vectors.

    EXAMPLES: 
    
    >> vecwt(2, 3, GF(2))
    [(1, 1, 0), (1, 0, 1), (0, 1, 1)]

    """
    L = []
    Kstar = K.list()[1:]
    fqw = cartesian_product([Kstar for i in range(w)])
    S = combinations(range(n), w) 
    for s in S:
        for i in fqw:
            vec = []
            counter = 0
            for j in range(n):
                if j in s:
                    vec.append(i[counter])
                    counter = counter + 1
                else:
                    vec.append(0)
            L.append(vector(K, vec))
    return L
    
def colwt(w, n, K):
    r"""
    This is an auxiliary function for GHW, hierarchy, RGHW and rhierarchy. 
    It computes all the vectors v of length n such that the last n - w coordinates
    are equal to 0 and Hamming weight between 1 and w (both included).

    OUTPUT:

    A list of such vectors.

    EXAMPLES:

    >> colwt(2, 3, GF(2))
    [(1, 0, 0), (0, 1, 0), (1, 1, 0)]
    
    """
    L = []
    for i in range(1, w + 1):
        L = L + vecwt(i, w, K)
    return [vector(K, list(i) + [0 for j in range(n - w)]) for i in L]
    
def standard(i, n, K):
    r"""
    Returns the standard vector (0,...,0,1,0,...,0) of length n with a 1 in the i-th
    coordinate and zeroes in the rest, over the finite field K.

    OUTPUT:

    The corresponding standard vector.

    EXAMPLES:

    >> standard(2, 3, GF(2))
    (0, 0, 1)
    
    """
    temp = [0 for j in range(i)] + [1] + [0 for j in range(i, n - 1)]
    return vector(K, temp)

def is_cyclic(C):
    r""""
    Returns True if the code C is cyclic, and False otherwise.

    OUTPUT:

    True if C is cyclic, False otherwise

    EXAMPLES:

    >> C = codes.BCHCode(GF(2), 15, 3)
    >> is_cyclic(C)
    True

    >> C = codes.BinaryReedMullerCode(2, 4)
    >> is_cyclic(C)
    False

    """
    G = C.generator_matrix()
    K = C.base_field()
    n = C.length()
    P = Permutation([n] + [i for i in range(1, n)])
    cyc = True
    for i in G:
        if vector(K, P.action(i)) not in C:
            cyc = False
            break
    return cyc

def bch_bound(C):
    r"""
    Returns the BCH bound for the minimum distance of C, 
    if C is a cyclic code.

    OUTPUT:

    The value of the bound.

    EXAMPLES:

    >> C = codes.BCHCode(GF(2), 15, 3)
    >> bch_bound(C)
    3

    """
    if not is_cyclic(C):
        raise Exception('The code is not cyclic')
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if k == n:
        return 1
    F = []
    R = PolynomialRing(K, 'x')
    x = R.gens()[0]
    for i in G:
        f = 0
        for j in range(n):
            f = f + i[j] * x**j
        F.append(f)
    g = gcd(F)
    CC = codes.CyclicCode(generator_pol = g, length = n)
    I = CC.defining_set()
    Ibis = I + I # We duplicate I to be able to detect consecutive elements containing
    # both 0 and n - 1
    bound = 1
    for i in range(len(I)):
        consecutive = 0
        for j in range(len(I)):
            if (I[i] + j in Ibis) % n: 
                consecutive = consecutive + 1
            else:
                break
        bound = max(bound, consecutive + 1)
    return bound

def information(G):
    r"""
    Computes a list of information sets inf_sets=[I1,...,Im] for the code generated
    by the matrix G that covers all the entries {1,...,n} (n is the length of C), 
    computes a list gen_mats of generator matrices for C that are systematic in each
    of the information sets in inf_sets, and computes the redundancies of the sets in
    inf_sets.

    OUTPUT:

    A list [inf_sets,gen_mats,redundancy] with three lists. The list inf_sets contains
    the information sets, gen_mats contains the corresponding generator matrices in 
    systematic form, and redundancy contains the values of the redundancies corresponding
    to the information sets in inf_sets.

    EXAMPLE:

    >> C = codes.BCHCode(GF(2), 15, 3)
    >> information(C.generator_matrix())
    [[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14]],
     [
    [1 0 0 0 0 0 0 0 0 0 0 1 1 0 0]  [0 0 0 0 0 0 0 1 1 0 0 1 0 0 0]
    [0 1 0 0 0 0 0 0 0 0 0 0 1 1 0]  [0 0 0 0 0 0 0 0 1 1 0 0 1 0 0]
    [0 0 1 0 0 0 0 0 0 0 0 0 0 1 1]  [0 0 0 0 0 0 0 0 0 1 1 0 0 1 0]
    [0 0 0 1 0 0 0 0 0 0 0 1 1 0 1]  [0 0 0 0 0 0 0 1 1 0 1 0 0 0 1]
    [0 0 0 0 1 0 0 0 0 0 0 1 0 1 0]  [1 0 0 0 0 0 0 1 0 1 0 0 0 0 0]
    [0 0 0 0 0 1 0 0 0 0 0 0 1 0 1]  [0 1 0 0 0 0 0 0 1 0 1 0 0 0 0]
    [0 0 0 0 0 0 1 0 0 0 0 1 1 1 0]  [0 0 1 0 0 0 0 1 1 1 0 0 0 0 0]
    [0 0 0 0 0 0 0 1 0 0 0 0 1 1 1]  [0 0 0 1 0 0 0 0 1 1 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 1 0 0 1 1 1 1]  [0 0 0 0 1 0 0 1 1 1 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 1 0 1 0 1 1]  [0 0 0 0 0 1 0 1 0 1 1 0 0 0 0]
    [0 0 0 0 0 0 0 0 0 0 1 1 0 0 1], [0 0 0 0 0 0 1 1 0 0 1 0 0 0 0]
    ],
     [0, 7]]

    """
    G = G.rref() 
    n = G.ncols()
    k = G.rank()
    K = G.base_ring()
    inf_sets = []
    gen_mats = []
    sup = list(matrix_supp(G))
    nn = len(sup)
    GG = copy(G)
    P = range(1, n + 1)
    comp = sup # This contains the complement of the coordinates already
    # covered by the information sets in inf_sets
    while len(flatten(inf_sets)) < nn:
        GG = matrix(K, [G.columns()[i] for i in comp]).transpose() # We restrict
        # to the columns of G not yet covered 

        # We compute an information set
        Itemp = list(GG.pivots())
        Itemp = [comp[i] for i in Itemp] 
        if len(Itemp) == k:
            inf_sets.append(Itemp)
        else: # If the corresponding matrix does not have rank k, we need 
            # to consider information sets with some redundancy 
            G0 = gen_mats[0].columns()
            cols = [G0[i] for i in range(n) if i in Itemp]
            for i in range(n):
                if i not in Itemp and matrix(K, cols + [G0[i]]).rank()>len(Itemp):
                    Itemp.append(i)
                    cols.append(G0[i])
                if len(Itemp) == k:
                    break
            inf_sets.append(sorted(Itemp))

        # We construct the corresponding generator matrix
        Itemp2 = [i + 1 for i in Itemp]
        # This permutation puts the indices in Itemp in the first positions
        P = Permutation(Itemp2 + [i for i in range(1, n + 1) if i not in Itemp2])
        # We put the columns associated to the indices of Itemp first, and 
        # compute the reduced row echelon form
        PG = matrix(K, P.action(G.columns())).transpose().rref()
        # We go back to the initial ordering of the columns
        PPG = matrix(K, P.inverse().action(PG.columns())).transpose()
        gen_mats.append(PPG)        
    
        I = flatten(inf_sets)
        comp = [i for i in sup if i not in I]
    # We compute the redundancies
    redundancy = [0] 
    for i in range(1, len(inf_sets)):
        previous = inf_sets[:i]
        previous = list(set([j for y in previous for j in y])) # ALl indices covered by
        # previous = inf_sets[:i]
        red = [y for y in inf_sets[i] if y in previous]
        redundancy.append(len(red))
    return [inf_sets, gen_mats, redundancy]

#################################################################################
# MAIN FUNCTIONS   
#################################################################################

def matrix_supp(M):
    r"""
    Returns the support of row space of the matrix M. This is equivalent to
    computing the union of the supports of the rows.

    OUTPUT:

    The support of the rowspace of M, as a set.

    EXAMPLES:

    >> M = matrix(GF(2), [[1,0,1,0],[1,1,0,0]])
    >> matrix_supp(M)
    {0,1,2}
    
    """
    sup = set()
    for i in M:
        sup.update(i.support())
    return sup

def subspaces(r, w, n, K):
    r"""
    Computes all the subspaces E contained in K^n with dim(E)=r and |supp(E)|=w.

    OUTPUT:

    A list with matrices in reduced row echelon form. The list of row spaces of
    these matrices gives all the aforementioned subspaces, without repetitions.

    EXAMPLES:

    >> subspaces(2, 3, 4, GF(2))
    [
    [1 0 1 0]  [1 0 0 0]  [1 0 1 0]  [1 1 0 0]  [1 0 0 1]  [1 0 0 0]
    [0 1 0 0], [0 1 1 0], [0 1 1 0], [0 0 1 0], [0 1 0 0], [0 1 0 1],
    
    [1 0 0 1]  [1 1 0 0]  [1 0 0 1]  [1 0 0 0]  [1 0 0 1]  [1 0 1 0]
    [0 1 0 1], [0 0 0 1], [0 0 1 0], [0 0 1 1], [0 0 1 1], [0 0 0 1],
    
    [0 1 0 1]  [0 1 0 0]  [0 1 0 1]  [0 1 1 0]
    [0 0 1 0], [0 0 1 1], [0 0 1 1], [0 0 0 1]
    ]
    
    """
    L = []
    for y in combinations(range(n), w):
        S = combinations(y[1:], r - 1) # We assume we have a pivot on the first
        # position, and we choose r - 1 more pivots
        for s in S:
            comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
            ss = [y[0]] + list(s) # Pivots
            comp = list(comp)
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
                MM = []
                for j in range(n):
                    if j in ss: # Pivots
                        MM.append(standard(ss.index(j), r, K))
                    elif j in comp: # Non pivots
                        MM.append(i[comp.index(j)])
                    else: # 0's outside of the support
                        MM.append(vector(K, [0 for zz in range(r)]))
                L.append(matrix(K, MM).transpose())
    return L

def num_subspaces(r, w, n, K):
    r"""
    Computes e_w^r, the number of subspaces E in K^n with dim(E)=r and
    |supp(E)|=w.

    OUTPUT:

    A non-negative integer.

    EXAMPLES:

    >> num_subspaces(2, 3, 4, GF(2))
    16

    """
    try:
        q = len(K.list())
    except:
        raise Exception('K has to be a finite field')
    s = 0
    for i in range(w - r + 1):
        s = s + (-1)**i * gaussian_binomial(w - i, r)(q=q) * binomial(w, i)
    return s * binomial(n, w)
    
def wei_duality(L, n=None): 
    r"""
    When provided with L, a list with the weight hierarchy of a code C, it
    computes the GHWs of the dual via Wei duality. The optional argument n 
    is the length of the code C. If none is given, we assume the length of C
    is equal to L[-1], i.e., C is non-degenerate.

    OUTPUT:

    A list with the weight hierarchy of the dual code.

    EXAMPLES:

    >> wei_duality([5,7,10,11,12,13,14,15])
    [7, 8, 10, 12, 13, 14, 15]

    """
    if n is None:
        n = L[-1]
    return [i for i in range(1, n + 1) if n - i + 1 not in L]
    
def GHW(C, r, L=None, verbose=False):  
    r"""
    Computes the rth GHW of the code C. The optional argument L is a list
    of the form [inf_sets, gen_mats, redundancy], where inf_sets is a list
    of information sets that covers {1,...,n}, gen_mats are the corresponding
    generator matrices of C that are systematic in those information sets, and
    redundancy is the list of the redundancies of the information sets. If none
    is given, we use L=information(G) for a generator matrix G of C (unless the 
    code is cyclic, in which case we use additional improvements). If we set 
    verbose=True, real time information is provided during the excution of the 
    algorithm (see the example).

    OUTPUT:

    The value of the rth GHW of C. 

    EXAMPLES:

    >> C = codes.BinaryReedMullerCode(1, 5)
    >> GHW(C, 2)
    24
    >> GHW(C, 2, verbose=True)
    Lower: 2 Upper: 28 Support: 2 Expected: 5
    Subspace with cardinality of support 24 found
    Lower: 14 Upper: 24 Support: 3 Expected: 4
    Lower: 19 Upper: 24 Support: 4 Expected: 4
    24
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception('Invalid value of r')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C)) 
        except:
            ghwlb = ghwlb
    ghwub = n - k + r 
    w = r
    while w <= k and ghwlb < ghwub: 
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < ghwub:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < ghwub:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= ghwub:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to get ghwlb >= ghwub at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
            
        # These are the only matrices that contribute       
        gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
        terminate = False
        rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
        for y in combinations(range(k), w): # All possible supports of weight w
            for mat in rrefs:
                cols = mat.columns()
                MM = []
                for t in range(k):
                    if t in y:
                        MM.append(cols[y.index(t)])
                    else:
                        MM.append([0 for z in range(r)])
                Mtemp = matrix(MM).transpose()
                for j in gen_reduced:
                    soptemp = len(matrix_supp(Mtemp * j))
                    if soptemp < ghwub:
                        ghwub = soptemp
                        if verbose:
                            print('Subspace with cardinality of support', soptemp, 'found')
                        if ghwub <= ghwlb:
                            terminate = True
                            break
                if terminate:
                    break
            if terminate:
                break
        # Lower bound calculations     
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        w = w + 1
    return ghwub

def hierarchy(C, L=None, verbose=False): 
    r"""
    Computes the weight hierarchy of C. The optional argument L is a list
    of the form [inf_sets, gen_mats, redundancy], where inf_sets is a list
    of information sets that covers {1,...,n}, gen_mats are the corresponding
    generator matrices of C that are systematic in those information sets, and
    redundancy is the list of the redundancies of the information sets. If none
    is given, we use L=information(G) for a generator matrix G of C (unless the 
    code is cyclic, in which case we use additional improvements). If we set 
    verbose=True, real time information is provided during the excution of the 
    algorithm (see the example).

    OUTPUT:

    A list with the weight hierarchy of C.

    EXAMPLES:
    >> C = codes.BinaryReedMullerCode(1, 5)
    >> hierarchy(C)
    [16, 24, 28, 30, 31, 32]
    >> hierarchy(C, verbose=True)
    r = 1
    Lower: 1 Upper: 27 Support: 1 Expected: 5
    Subspace with cardinality of support 16 found
    Lower: 9 Upper: 16 Support: 2 Expected: 3
    Lower: 12 Upper: 16 Support: 3 Expected: 3
    Current hierarchy: [16]
    r = 2
    Lower: 17 Upper: 28 Support: 2 Expected: 5
    Subspace with cardinality of support 24 found
    Lower: 17 Upper: 24 Support: 3 Expected: 4
    Lower: 19 Upper: 24 Support: 4 Expected: 4
    Current hierarchy: [16, 24]
    r = 3
    Lower: 25 Upper: 29 Support: 3 Expected: 5
    Subspace with cardinality of support 28 found
    Lower: 25 Upper: 28 Support: 4 Expected: 5
    Lower: 25 Upper: 28 Support: 5 Expected: 5
    Current hierarchy: [16, 24, 28]
    r = 4
    Lower: 29 Upper: 30 Support: 4 Expected: 5
    Lower: 29 Upper: 30 Support: 5 Expected: 5
    Current hierarchy: [16, 24, 28, 30]
    r = 5
    Current hierarchy: [16, 24, 28, 30, 31]
    r = 6
    Current hierarchy: [16, 24, 28, 30, 31, 32]
    [16, 24, 28, 30, 31, 32]

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    hierarchyghw = []
    for r in range(1, k + 1):
        if verbose:
            print('r =', r)
        if r == 1:
            ghwlb = r
            if cyc:
                try:
                    ghwlb = max(ghwlb, bch_bound(C)) 
                except:
                    ghwlb = ghwlb
        else:
            ghwlb = hierarchyghw[-1] + 1 
        ghwub = n - k + r 
        w = r
        while w <= k and ghwlb < ghwub: 
            terminate = False
            # Computation of w0, the expected w to finish
            rm = []
            if cyc:
                w0 = w
                while ceil((w0 + 1) * n / k) < ghwub:
                    w0 = w0 + 1
            else:
                w0 = w - 1
                ghwlbtemp = 0
                while ghwlbtemp < ghwub:
                    ghwlbtemp = 0
                    w0 = w0 + 1
                    for j in range(len(gen)):
                        ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                        if ghwlbtemp >= ghwub:
                            if w0 == w:
                                # We store in rm the indices of the matrices that are not necessary
                                # to get ghwlb >= ghwub at the end of this iteration (if any)
                                rm = rm + srange(j + 1, len(gen))
                            break
            if verbose:
                print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
                
            # These are the only matrices that contribute       
            gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
            terminate = False
            rrefs = subspaces(r, w, w, K)  # All reduced row echelon forms in the first w columns
            for y in combinations(range(k), w): # All possible supports of weight w
                for mat in rrefs:
                    cols = mat.columns()
                    MM = []
                    for t in range(k):
                        if t in y:
                            MM.append(cols[y.index(t)])
                        else:
                            MM.append([0 for z in range(r)])
                    Mtemp = matrix(MM).transpose()
                    for j in gen_reduced:
                        soptemp = len(matrix_supp(Mtemp * j))
                        if soptemp < ghwub:
                            ghwub = soptemp
                            if verbose:
                                print('Subspace with cardinality of support', soptemp, 'found')
                            if ghwub <= ghwlb:
                                terminate = True
                                break                    
                    if terminate:
                        break
                if terminate:
                    break
            # Lower bound calculations        
            ghwlbtemp = 0
            for j in range(len(gen_reduced)):
                if cyc:
                    ghwlbtemp = ceil((w + 1) * n / k)
                else:
                    ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
            ghwlb = max(ghwlbtemp, ghwlb)
            w = w + 1
        hierarchyghw.append(ghwub)
        if verbose:
            print('Current hierarchy:', hierarchyghw)
    return hierarchyghw

def RGHW(C, C2, r, L=None, verbose=False):
    r"""
    Computes the rth RGHW of the code C with respect to the code C2. 
    The optional argument L is a list of the form [inf_sets, gen_mats, redundancy], 
    where inf_sets is a list of information sets for C that covers {1,...,n}, 
    gen_mats are the corresponding generator matrices of C that are systematic
    in those information sets, and redundancy is the list of the redundancies 
    of the information sets. If none is given, we use L=information(G) for a 
    generator matrix G of C (unless the code is cyclic, in which case we use 
    additional improvements). If we set verbose=True, real time information is 
    provided during the excution of the algorithm (see the example).

    OUTPUT:

    The value of the rth RGHW of C with respect to C2. 

    EXAMPLES:

    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> RGHW(C, C2, 2)
    5
    >> RGHW(C, C2, 2, verbose=True)
    Lower: 2 Upper: 6 Support: 2 Expected: 2
    Subspace with cardinality of support 5 found
    5
    >> # Note that RGHW(C, C2, 2) can be higher than GHW(C, 2):
    >> GHW(C, 2)
    4
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if r not in range(1, k - k2 + 1):
        raise Exception('Invalid value of r')
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C)) 
        except:
            ghwlb = ghwlb
    ghwub = n - k + r
    w = r
    while w <= k and ghwlb < ghwub: 
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < ghwub:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < ghwub:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp+max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= ghwub:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to get ghwlb >= ghwub at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
            
        # These are the only matrices that contribute       
        gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
        terminate = False
        rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
        for y in combinations(range(k), w): # All possible supports of weight w
            for mat in rrefs:
                cols = mat.columns()
                MM = []
                for t in range(k):
                    if t in y:
                        MM.append(cols[y.index(t)])
                    else:
                        MM.append([0 for z in range(r)])
                Mtemp = matrix(MM).transpose()
                for j in gen_reduced:
                    Mtempj = Mtemp * j
                    soptemp = len(matrix_supp(Mtempj))
                    if soptemp < ghwub:
                        interd = r - (H2 * Mtempj.transpose()).rank()
                        if interd == 0:
                            if verbose:
                                print('Subspace with cardinality of support', soptemp, 'found')
                            ghwub = soptemp
                            if ghwub <= ghwlb:
                                terminate = True
                                break
                if terminate:
                    break
            if terminate:
                break
        # Lower bound calculations     
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        w = w + 1
    return ghwub

def rhierarchy(C, C2, L=None, verbose=False):  
    r"""
    Computes the relative weight hierarchy of the code C with respect to the code C2. 
    The optional argument L is a list of the form [inf_sets, gen_mats, redundancy], 
    where inf_sets is a list of information sets for C that covers {1,...,n}, 
    gen_mats are the corresponding generator matrices of C that are systematic
    in those information sets, and redundancy is the list of the redundancies 
    of the information sets. If none is given, we use L=information(G) for a 
    generator matrix G of C (unless the code is cyclic, in which case we use 
    additional improvements). If we set verbose=True, real time information is 
    provided during the excution of the algorithm (see the example).

    OUTPUT:

    A list with the relative weight hierarchy of C with respect to C2.

    EXAMPLES:

    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> rhierarchy(C, C2)
    [3, 5]
    >> rhierarchy(C, C2, verbose=True)
    r = 1
    Lower: 1 Upper: 6 Support: 1 Expected: 2
    Subspace with cardinality of support 3 found
    Current hierarchy: [3]
    r = 2
    Lower: 4 Upper: 6 Support: 2 Expected: 2
    Subspace with cardinality of support 5 found
    Current hierarchy: [3, 5]
    [3, 5]
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    hierarchyrghw = []
    for r in range(1, k - k2 + 1):
        if verbose:
            print('r =', r)
        if r == 1:
            ghwlb = r
            if cyc:
                try:
                    ghwlb = max(ghwlb, bch_bound(C)) 
                except:
                    ghwlb = ghwlb
        else:
            ghwlb = hierarchyrghw[-1] + 1 
        ghwub = n - k + r
        w = r
        while w <= k and ghwlb < ghwub: 
            # Computation of w0, the expected w to finish
            rm = []
            if cyc:
                w0 = w
                while ceil((w0 + 1) * n / k) < ghwub:
                    w0 = w0 + 1
            else:
                w0 = w - 1
                ghwlbtemp = 0
                while ghwlbtemp < ghwub:
                    ghwlbtemp = 0
                    w0 = w0 + 1
                    for j in range(len(gen)):
                        ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                        if ghwlbtemp >= ghwub:
                            if w0 == w:
                                # We store in rm the indices of the matrices that are not necessary
                                # to get ghwlb >= ghwub at the end of this iteration (if any)
                                rm = rm + srange(j + 1, len(gen))
                            break
            if verbose:
                print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
                
            # These are the only matrices that contribute       
            gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
            terminate = False
            rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
            for y in combinations(range(k), w): # All possible supports of weight w
                for mat in rrefs:
                    cols = mat.columns()
                    MM = []
                    for t in range(k):
                        if t in y:
                            MM.append(cols[y.index(t)])
                        else:
                            MM.append([0 for z in range(r)])
                    Mtemp = matrix(MM).transpose()
                    for j in gen_reduced:
                        Mtempj = Mtemp * j
                        soptemp = len(matrix_supp(Mtempj))
                        if soptemp < ghwub:
                            interd = r - (H2 * Mtempj.transpose()).rank()
                            if interd == 0:
                                if verbose:
                                    print('Subspace with cardinality of support', soptemp, 'found')
                                ghwub = soptemp
                                if ghwub <= ghwlb:
                                    terminate = True
                                    break
                    if terminate:
                        break
                if terminate:
                    break
            # Lower bound calculations        
            ghwlbtemp = 0
            for j in range(len(gen_reduced)):
                if cyc:
                    ghwlbtemp = ceil((w + 1) * n / k)
                else:
                    ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
            ghwlb = max(ghwlbtemp, ghwlb)
            w = w + 1
        hierarchyrghw.append(ghwub)
        if verbose:
            print('Current hierarchy:', hierarchyrghw)
    return hierarchyrghw

def higher_spectrum(C, verbose=False, subspace_list=False, R=None): 
    r"""
    Computes the higher weight spectrum of the code C. If we set verbose=True, 
    real time information is provided during the excution of the algorithm
    (see the example). If we set subspace_list=True, the function will also
    return a dictionary 'subsp', whose keys are equal to the corresponding r,
    and whose values are dictionaries, whose keys are equal to each possible 
    cardinality of the support w, and whose values are the lists of all subcodes 
    of dimension r and cardinality of support w (given as matrices, the corresponding
    subspaces are the row spaces). In other words, subsp[r][w] is a list of the
    subspaces D in C with dim(D)=r, |supp(D)|=w.
    The optional parameter R allows the user to only compute the higher weight 
    spectrum for r in R. If none is given, we compute it for 1 <= r <= dim(C). 

    OUTPUT:

    If subspace_list=False, a dictionary 'spec', whose keys are given by r, and
    whose values are dictionaries, whose keys are the possible weights w, and
    whose values are the number of subcodes D of C with dim(D)=r and |supp(E)|=w.
    If subspace_list=False, a list [spec, subsp], where 'spec' and 'subsp'
    are the dictionaries mentioned above. 

    EXAMPLES:

    >> C = codes.GeneralizedReedSolomonCode(GF(7).list(), 4)
    >> higher_spectrum(C)
    {1: {4: 35, 5: 63, 6: 168, 7: 134},
     2: {5: 21, 6: 357, 7: 2472},
     3: {6: 7, 7: 393},
     4: {7: 1}}
    >> higher_spectrum(C, verbose=True) 
    r = 1
    Current spectrum: {4: 35, 5: 63, 6: 168, 7: 134}
    r = 2
    Current spectrum: {5: 21, 6: 357, 7: 2472}
    r = 3
    Current spectrum: {6: 7, 7: 393}
    r = 4
    Current spectrum: {7: 1}
    {1: {4: 35, 5: 63, 6: 168, 7: 134},
     2: {5: 21, 6: 357, 7: 2472},
     3: {6: 7, 7: 393},
     4: {7: 1}}
    >> higher_spectrum(C, R=[2,3]) 
    {2: {5: 21, 6: 357, 7: 2472}, 3: {6: 7, 7: 393}}
    >> higher_spectrum(C, R=[4], subspace_list=True) 
    [{4: {7: 1}},
     {4: {7: [
    [1 0 0 0 6 3 4]
    [0 1 0 0 4 1 1]
    [0 0 1 0 1 1 4]
    [0 0 0 1 4 3 6]
    ]}}]
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    spec = {} 
    if subspace_list:
        subsp = {} 
    if R is None:
        R = range(k + 1)     
    for r in R:
        if verbose:
            print('r =', r)
        if r == 0:
            spectemp = {0: 1}
            spec.update({0: spectemp})
            if subspace_list:
                subsptemp = {0: [matrix(K, [0 for i in range(n)])]}
                subsp = {0: subsptemp}
            if verbose:
                print('Spectrum:', spec[r]) 
            continue 
        spectemp = {}
        if subspace_list:
            subsptemp = {}
        for w in range(r, k + 1):
            rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
            for y in combinations(range(k), w): # All possible supports of weight w
                for mat in rrefs:
                    cols = mat.columns()
                    MM = []
                    for t in range(k):
                        if t in y:
                            MM.append(cols[y.index(t)])
                        else:
                            MM.append([0 for z in range(r)])
                    Mtemp = matrix(MM).transpose() * G
                    soptemp = len(matrix_supp(Mtemp))
                    try:
                        spectemp.update({soptemp: spectemp[soptemp] + 1})
                    except:
                        spectemp.update({soptemp: 1})
                    if subspace_list:
                        try:
                            subsptemp.update({soptemp: subsptemp[soptemp] + [Mtemp]})
                        except:
                            subsptemp.update({soptemp: [Mtemp]})
            if verbose:
                spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} # Order by key
                print('Current spectrum:', spectemp, 'w =', w, end='\r') 
        spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} 
        spec.update({r: spectemp})
        if subspace_list:
            subsptemp = {k: v for k, v in sorted(subsptemp.items(), key=lambda item: item[0])} 
            subsp.update({r: subsptemp})
        if verbose:
            print('Spectrum:', spec[r], ' '*100) 
    if subspace_list:
        return [spec, subsp]
    else:
        return spec
    
def rhigher_spectrum(C, C2, verbose=False, subspace_list=False, R=None): 
    r"""
    Computes the relative higher weight spectrum of the code C. If we set 
    verbose=True, real time information is provided during the excution of 
    the algorithm (see the example). If we set subspace_list=True, the function 
    will also return a dictionary 'subsp', whose keys are equal to the 
    corresponding r, and whose values are dictionaries, whose keys are equal to 
    each possible cardinality of the support w, and whose values are the lists of 
    all subcodes of dimension r and cardinality of support w (given as matrices, 
    the corresponding subspaces are the row spaces). In other words, subsp[r][w] 
    is a list of the subspaces D in C with dim(D)=r, |supp(D)|=w, and with trivial
    intersection with C2.
    The optional parameter R allows the user to only compute the relative higher weight 
    spectrum for r in R. If none is given, we compute it for 1 <= r <= dim(C)-dim(C2). 
    
    OUTPUT:

    If subspace_list=False, a dictionary 'spec', whose keys are given by r, and
    whose values are dictionaries, whose keys are the possible weights w, and
    whose values are the number of subcodes D of C with dim(D)=r and |supp(E)|=w.
    If subspace_list=False, a list [spec, subsp], where 'spec' and 'subsp'
    are the dictionaries mentioned above. 

    EXAMPLES:
    
    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> rhigher_spectrum(C, C2)
    {0: {0: 1}, 1: {3: 4, 4: 1, 6: 1}, 2: {5: 2, 6: 2}}
    >> rhigher_spectrum(C, C2, verbose=True)
    r = 0
    Spectrum: {0: 1}
    r = 1
    Spectrum: {3: 4, 4: 1, 6: 1}                                                                                                     
    r = 2
    Spectrum: {5: 2, 6: 2}                                                                                                     
    {0: {0: 1}, 1: {3: 4, 4: 1, 6: 1}, 2: {5: 2, 6: 2}}
    >> rhigher_spectrum(C, C2, R=[1])
    {1: {3: 4, 4: 1, 6: 1}}
    >> rhigher_spectrum(C, C2, R=[2], subspace_list=True)
    [{2: {5: 2, 6: 2}},
     {2: {5: [
    [1 0 0 0 1 1]  [1 0 1 0 0 1]
    [0 1 0 1 1 0], [0 1 1 1 0 0]
    ],
       6: [
    [1 0 1 0 0 1]  [1 0 0 0 1 1]
    [0 1 0 1 1 0], [0 1 1 1 0 0]
    ]}}]

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    spec = {} 
    if subspace_list:
        subsp = {} 
    if R is None:
        R = range(k - k2 + 1) 
    elif not set(R).issubset(set(range(k - k2 + 1))):
        raise Exception('Invalid range of values for R')
        
    for r in R:
        if verbose:
            print('r =', r)
        if r == 0:
            spectemp = {0: 1}
            spec.update({0: spectemp})
            if subspace_list:
                subsptemp = {0: [matrix(K, [0 for i in range(n)])]}
                subsp = {0: subsptemp}
            if verbose:
                print('Spectrum:', spec[r]) 
            continue  
        spectemp = {}
        if subspace_list:
            subsptemp = {}
        for w in range(r, k + 1):
            rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
            for y in combinations(range(k), w): # All possible supports of weight w
                for mat in rrefs:
                    cols = mat.columns()
                    MM = []
                    for t in range(k):
                        if t in y:
                            MM.append(cols[y.index(t)])
                        else:
                            MM.append([0 for z in range(r)])
                    Mtemp = matrix(MM).transpose() * G
                    soptemp = len(matrix_supp(Mtemp))
                    interd = r - (H2 * Mtemp.transpose()).rank()
                    if interd == 0:
                        try:
                            spectemp.update({soptemp: spectemp[soptemp] + 1})
                        except:
                            spectemp.update({soptemp: 1})
                        if subspace_list:
                            try:
                                subsptemp.update({soptemp: subsptemp[soptemp] + [Mtemp]})
                            except:
                                subsptemp.update({soptemp: [Mtemp]})
            if verbose:
                spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} # Order by key
                print('Current spectrum:', spectemp, 'w =', w, end='\r') 
        spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} 
        spec.update({r: spectemp})
        if subspace_list:
            subsptemp = {k: v for k, v in sorted(subsptemp.items(), key=lambda item: item[0])} 
            subsp.update({r: subsptemp})
        if verbose:
            print('Spectrum:', spec[r], ' '*100) 
    if subspace_list:
        return [spec, subsp]
    else:
        return spec
    

# Low memory functions

def GHW_low_mem(C, r, L=None, verbose=False):  
    r"""
    Computes the rth GHW of the code C. The optional argument L is a list
    of the form [inf_sets, gen_mats, redundancy], where inf_sets is a list
    of information sets that covers {1,...,n}, gen_mats are the corresponding
    generator matrices of C that are systematic in those information sets, and
    redundancy is the list of the redundancies of the information sets. If none
    is given, we use L=information(G) for a generator matrix G of C (unless the 
    code is cyclic, in which case we use additional improvements). If we set 
    verbose=True, real time information is provided during the excution of the 
    algorithm (see the example). This is a version of GHW that requires less
    memory, at the expense of speed in some cases.

    OUTPUT:

    The value of the rth GHW of C. 

    EXAMPLES:

    >> C = codes.BinaryReedMullerCode(1, 5)
    >> GHW_low_mem(C, 2)
    24
    >> GHW_low_mem(C, 2, verbose=True)
    Lower: 2 Upper: 28 Support: 2 Expected: 5
    Subspace with cardinality of support 24 found
    Lower: 14 Upper: 24 Support: 3 Expected: 4
    Lower: 19 Upper: 24 Support: 4 Expected: 4
    24
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception('Invalid value of r')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C)) 
        except:
            ghwlb = ghwlb
    ghwub = n - k + r 
    w = r
    while w <= k and ghwlb < ghwub: 
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < ghwub:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < ghwub:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= ghwub:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to get ghwlb >= ghwub at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
            
        # These are the only matrices that contribute       
        gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
        terminate = False
        y = range(w) # We start with support {1,...,w}
        for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
            # position, and we choose r - 1 more pivots
            comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
            ss = [y[0]] + list(s) # Pivots
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
                    
                    Mtemp = matrix(K, MM).transpose()

                    for j in gen_reduced:
                        soptemp = len(matrix_supp(Mtemp * j))
                        if soptemp < ghwub:
                            ghwub = soptemp
                            if verbose:
                                print('Subspace with cardinality of support', soptemp, 'found')
                            if ghwub <= ghwlb:
                                terminate = True
                                break
                    if terminate:
                        break
                if terminate:
                    break
            if terminate:
                break
        # Lower bound calculations     
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        w = w + 1
    return ghwub

def hierarchy_low_mem(C, L=None, verbose=False): 
    r"""
    Computes the weight hierarchy of C. The optional argument L is a list
    of the form [inf_sets, gen_mats, redundancy], where inf_sets is a list
    of information sets that covers {1,...,n}, gen_mats are the corresponding
    generator matrices of C that are systematic in those information sets, and
    redundancy is the list of the redundancies of the information sets. If none
    is given, we use L=information(G) for a generator matrix G of C (unless the 
    code is cyclic, in which case we use additional improvements). If we set 
    verbose=True, real time information is provided during the excution of the 
    algorithm (see the example). This is a version of hierarchy that requires 
    less memory, at the expense of speed in some cases.

    OUTPUT:

    A list with the weight hierarchy of C.

    EXAMPLES:
    >> C = codes.BinaryReedMullerCode(1, 5)
    >> hierarchy_low_mem(C)
    [16, 24, 28, 30, 31, 32]
    >> hierarchy_low_mem(C, verbose=True)
    r = 1
    Lower: 1 Upper: 27 Support: 1 Expected: 5
    Subspace with cardinality of support 16 found
    Lower: 9 Upper: 16 Support: 2 Expected: 3
    Lower: 12 Upper: 16 Support: 3 Expected: 3
    Current hierarchy: [16]
    r = 2
    Lower: 17 Upper: 28 Support: 2 Expected: 5
    Subspace with cardinality of support 24 found
    Lower: 17 Upper: 24 Support: 3 Expected: 4
    Lower: 19 Upper: 24 Support: 4 Expected: 4
    Current hierarchy: [16, 24]
    r = 3
    Lower: 25 Upper: 29 Support: 3 Expected: 5
    Subspace with cardinality of support 28 found
    Lower: 25 Upper: 28 Support: 4 Expected: 5
    Lower: 25 Upper: 28 Support: 5 Expected: 5
    Current hierarchy: [16, 24, 28]
    r = 4
    Lower: 29 Upper: 30 Support: 4 Expected: 5
    Lower: 29 Upper: 30 Support: 5 Expected: 5
    Current hierarchy: [16, 24, 28, 30]
    r = 5
    Current hierarchy: [16, 24, 28, 30, 31]
    r = 6
    Current hierarchy: [16, 24, 28, 30, 31, 32]
    [16, 24, 28, 30, 31, 32]

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    hierarchyghw = []
    for r in range(1, k + 1):
        if verbose:
            print('r =', r)
        if r == 1:
            ghwlb = r
            if cyc:
                try:
                    ghwlb = max(ghwlb, bch_bound(C)) 
                except:
                    ghwlb = ghwlb
        else:
            ghwlb = hierarchyghw[-1] + 1 
        ghwub = n - k + r 
        w = r
        while w <= k and ghwlb < ghwub: 
            terminate = False
            # Computation of w0, the expected w to finish
            rm = []
            if cyc:
                w0 = w
                while ceil((w0 + 1) * n / k) < ghwub:
                    w0 = w0 + 1
            else:
                w0 = w - 1
                ghwlbtemp = 0
                while ghwlbtemp < ghwub:
                    ghwlbtemp = 0
                    w0 = w0 + 1
                    for j in range(len(gen)):
                        ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                        if ghwlbtemp >= ghwub:
                            if w0 == w:
                                # We store in rm the indices of the matrices that are not necessary
                                # to get ghwlb >= ghwub at the end of this iteration (if any)
                                rm = rm + srange(j + 1, len(gen))
                            break
            if verbose:
                print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
                
            # These are the only matrices that contribute       
            gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
            terminate = False
            y = range(w) # We start with support {1,...,w}
            for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
                # position, and we choose r - 1 more pivots
                comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
                ss = [y[0]] + list(s) # Pivots
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
                        
                        Mtemp = matrix(K, MM).transpose()

                        for j in gen_reduced:
                            soptemp = len(matrix_supp(Mtemp * j))
                            if soptemp < ghwub:
                                ghwub = soptemp
                                if verbose:
                                    print('Subspace with cardinality of support', soptemp, 'found')
                                if ghwub <= ghwlb:
                                    terminate = True
                                    break
                        if terminate:
                            break
                    if terminate:
                        break
                if terminate:
                    break

            # Lower bound calculations        
            ghwlbtemp = 0
            for j in range(len(gen_reduced)):
                if cyc:
                    ghwlbtemp = ceil((w + 1) * n / k)
                else:
                    ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
            ghwlb = max(ghwlbtemp, ghwlb)
            w = w + 1
        hierarchyghw.append(ghwub)
        if verbose:
            print('Current hierarchy:', hierarchyghw)
    return hierarchyghw

def RGHW_low_mem(C, C2, r, L=None, verbose=False):
    r"""
    Computes the rth RGHW of the code C with respect to the code C2. 
    The optional argument L is a list of the form [inf_sets, gen_mats, redundancy], 
    where inf_sets is a list of information sets for C that covers {1,...,n}, 
    gen_mats are the corresponding generator matrices of C that are systematic
    in those information sets, and redundancy is the list of the redundancies 
    of the information sets. If none is given, we use L=information(G) for a 
    generator matrix G of C (unless the code is cyclic, in which case we use 
    additional improvements). If we set verbose=True, real time information is 
    provided during the excution of the algorithm (see the example). This is a 
    version of RGHW that requires less memory, at the expense of speed in some 
    cases.

    OUTPUT:

    The value of the rth RGHW of C with respect to C2. 

    EXAMPLES:

    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> RGHW_low_mem(C, C2, 2)
    5
    >> RGHW_low_mem(C, C2, 2, verbose=True)
    Lower: 2 Upper: 6 Support: 2 Expected: 2
    Subspace with cardinality of support 5 found
    5
    >> # Note that RGHW_low_mem(C, C2, 2) can be higher than GHW(C, 2):
    >> GHW(C, 2)
    4
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if r not in range(1, k - k2 + 1):
        raise Exception('Invalid value of r')
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]] 
        else:
            L = information(G)
    [inf, gen, red] = L
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C)) 
        except:
            ghwlb = ghwlb
    ghwub = n - k + r
    w = r
    while w <= k and ghwlb < ghwub: 
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < ghwub:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < ghwub:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp+max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= ghwub:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to get ghwlb >= ghwub at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
            
        # These are the only matrices that contribute 
        gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
        terminate = False
        y = range(w) # We start with support {1,...,w}
        for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
            # position, and we choose r - 1 more pivots
            comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
            ss = [y[0]] + list(s) # Pivots
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
                    
                    Mtemp = matrix(K, MM).transpose()

                    for j in gen_reduced:
                        Mtempj = Mtemp * j
                        soptemp = len(matrix_supp(Mtempj))
                        if soptemp < ghwub:
                            interd = r - (H2 * Mtempj.transpose()).rank()
                            if interd == 0:
                                if verbose:
                                    print('Subspace with cardinality of support', soptemp, 'found')
                                ghwub = soptemp
                                if ghwub <= ghwlb:
                                    terminate = True
                                    break
                    if terminate:
                        break
                if terminate:
                    break
            if terminate:
                break

        # Lower bound calculations     
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        w = w + 1
    return ghwub

def rhierarchy_low_mem(C, C2, L=None, verbose=False):  
    r"""
    Computes the relative weight hierarchy of the code C with respect to the code C2. 
    The optional argument L is a list of the form [inf_sets, gen_mats, redundancy], 
    where inf_sets is a list of information sets for C that covers {1,...,n}, 
    gen_mats are the corresponding generator matrices of C that are systematic
    in those information sets, and redundancy is the list of the redundancies 
    of the information sets. If none is given, we use L=information(G) for a 
    generator matrix G of C (unless the code is cyclic, in which case we use 
    additional improvements). If we set verbose=True, real time information is 
    provided during the excution of the algorithm (see the example). This is a 
    version of rhierarchy that requires less memory, at the expense of speed in 
    some cases.

    OUTPUT:

    A list with the relative weight hierarchy of C with respect to C2.

    EXAMPLES:

    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> rhierarchy_low_mem(C, C2)
    [3, 5]
    >> rhierarchy_low_mem(C, C2, verbose=True)
    r = 1
    Lower: 1 Upper: 6 Support: 1 Expected: 2
    Subspace with cardinality of support 3 found
    Current hierarchy: [3]
    r = 2
    Lower: 4 Upper: 6 Support: 2 Expected: 2
    Subspace with cardinality of support 5 found
    Current hierarchy: [3, 5]
    [3, 5]
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    # Only cyclic codes with non-repeted roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    hierarchyrghw = []
    for r in range(1, k - k2 + 1):
        if verbose:
            print('r =', r)
        if r == 1:
            ghwlb = r
            if cyc:
                try:
                    ghwlb = max(ghwlb, bch_bound(C)) 
                except:
                    ghwlb = ghwlb
        else:
            ghwlb = hierarchyrghw[-1] + 1 
        ghwub = n - k + r
        w = r
        while w <= k and ghwlb < ghwub: 
            # Computation of w0, the expected w to finish
            rm = []
            if cyc:
                w0 = w
                while ceil((w0 + 1) * n / k) < ghwub:
                    w0 = w0 + 1
            else:
                w0 = w - 1
                ghwlbtemp = 0
                while ghwlbtemp < ghwub:
                    ghwlbtemp = 0
                    w0 = w0 + 1
                    for j in range(len(gen)):
                        ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                        if ghwlbtemp >= ghwub:
                            if w0 == w:
                                # We store in rm the indices of the matrices that are not necessary
                                # to get ghwlb >= ghwub at the end of this iteration (if any)
                                rm = rm + srange(j + 1, len(gen))
                            break
            if verbose:
                print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)
                
            # These are the only matrices that contribute   
            gen_reduced = [gen[j] for j in range(len(gen)) if red[j] <= w and j not in rm] 
            terminate = False
            y = range(w) # We start with support {1,...,w}
            for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
                # position, and we choose r - 1 more pivots
                comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
                ss = [y[0]] + list(s) # Pivots
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
                        
                        Mtemp = matrix(K, MM).transpose()

                        for j in gen_reduced:
                            Mtempj = Mtemp * j
                            soptemp = len(matrix_supp(Mtempj))
                            if soptemp < ghwub:
                                interd = r - (H2 * Mtempj.transpose()).rank()
                                if interd == 0:
                                    if verbose:
                                        print('Subspace with cardinality of support', soptemp, 'found')
                                    ghwub = soptemp
                                    if ghwub <= ghwlb:
                                        terminate = True
                                        break
                        if terminate:
                            break
                    if terminate:
                        break
                if terminate:
                    break

            # Lower bound calculations        
            ghwlbtemp = 0
            for j in range(len(gen_reduced)):
                if cyc:
                    ghwlbtemp = ceil((w + 1) * n / k)
                else:
                    ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
            ghwlb = max(ghwlbtemp, ghwlb)
            w = w + 1
        hierarchyrghw.append(ghwub)
        if verbose:
            print('Current hierarchy:', hierarchyrghw)
    return hierarchyrghw

def higher_spectrum_low_mem(C, verbose=False, subspace_list=False, R=None): 
    r"""
    Computes the higher weight spectrum of the code C. If we set verbose=True, 
    real time information is provided during the excution of the algorithm
    (see the example). If we set subspace_list=True, the function will also
    return a dictionary 'subsp', whose keys are equal to the corresponding r,
    and whose values are dictionaries, whose keys are equal to each possible 
    cardinality of the support w, and whose values are the lists of all subcodes 
    of dimension r and cardinality of support w (given as matrices, the corresponding
    subspaces are the row spaces). In other words, subsp[r][w] is a list of the
    subspaces D in C with dim(D)=r, |supp(D)|=w.
    The optional parameter R allows the user to only compute the higher weight 
    spectrum for r in R. If none is given, we compute it for 1 <= r <= dim(C). 
    This is a version of higher_spectrum that requires less memory, at the expense 
    of speed in some cases.

    OUTPUT:

    If subspace_list=False, a dictionary 'spec', whose keys are given by r, and
    whose values are dictionaries, whose keys are the possible weights w, and
    whose values are the number of subcodes D of C with dim(D)=r and |supp(E)|=w.
    If subspace_list=False, a list [spec, subsp], where 'spec' and 'subsp'
    are the dictionaries mentioned above. 

    EXAMPLES:

    >> C = codes.GeneralizedReedSolomonCode(GF(7).list(), 4)
    >> higher_spectrum_low_mem(C)
    {1: {4: 35, 5: 63, 6: 168, 7: 134},
     2: {5: 21, 6: 357, 7: 2472},
     3: {6: 7, 7: 393},
     4: {7: 1}}
    >> higher_spectrum_low_mem(C, verbose=True) 
    r = 1
    Current spectrum: {4: 35, 5: 63, 6: 168, 7: 134}
    r = 2
    Current spectrum: {5: 21, 6: 357, 7: 2472}
    r = 3
    Current spectrum: {6: 7, 7: 393}
    r = 4
    Current spectrum: {7: 1}
    {1: {4: 35, 5: 63, 6: 168, 7: 134},
     2: {5: 21, 6: 357, 7: 2472},
     3: {6: 7, 7: 393},
     4: {7: 1}}
    >> higher_spectrum_low_mem(C, R=[2,3]) 
    {2: {5: 21, 6: 357, 7: 2472}, 3: {6: 7, 7: 393}}
    >> higher_spectrum_low_mem(C, R=[4], subspace_list=True) 
    [{4: {7: 1}},
     {4: {7: [
    [1 0 0 0 6 3 4]
    [0 1 0 0 4 1 1]
    [0 0 1 0 1 1 4]
    [0 0 0 1 4 3 6]
    ]}}]
    
    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    spec = {} 
    if subspace_list:
        subsp = {} 
    if R is None:
        R = range(k + 1)     
    for r in R:
        if verbose:
            print('r =', r)
        if r == 0:
            spectemp = {0: 1}
            spec.update({0: spectemp})
            if subspace_list:
                subsptemp = {0: [matrix(K, [0 for i in range(n)])]}
                subsp = {0: subsptemp}
            if verbose:
                print('Spectrum:', spec[r]) 
            continue 
        spectemp = {}
        if subspace_list:
            subsptemp = {}
        for w in range(r, k + 1):
            y = range(w) # We start with support {1,...,w}
            for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
                # position, and we choose r - 1 more pivots
                comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
                ss = [y[0]] + list(s) # Pivots
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
                        
                        Mtemp = matrix(K, MM).transpose() * G
                        soptemp = len(matrix_supp(Mtemp))
                        try:
                            spectemp.update({soptemp: spectemp[soptemp] + 1})
                        except:
                            spectemp.update({soptemp: 1})
                        if subspace_list:
                            try:
                                subsptemp.update({soptemp: subsptemp[soptemp] + [Mtemp]})
                            except:
                                subsptemp.update({soptemp: [Mtemp]})
            if verbose:
                spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} # Order by key
                print('Current spectrum:', spectemp, 'w =', w, end='\r') 
        spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} 
        spec.update({r: spectemp})
        if subspace_list:
            subsptemp = {k: v for k, v in sorted(subsptemp.items(), key=lambda item: item[0])} 
            subsp.update({r: subsptemp})
        if verbose:
            print('Spectrum:',spec[r], ' '*100) 
    if subspace_list:
        return [spec, subsp]
    else:
        return spec
    
def rhigher_spectrum_low_mem(C, C2, verbose=False, subspace_list=False, R=None): 
    r"""
    Computes the relative higher weight spectrum of the code C. If we set 
    verbose=True, real time information is provided during the excution of 
    the algorithm (see the example). If we set subspace_list=True, the function 
    will also return a dictionary 'subsp', whose keys are equal to the 
    corresponding r, and whose values are dictionaries, whose keys are equal to 
    each possible cardinality of the support w, and whose values are the lists of 
    all subcodes of dimension r and cardinality of support w (given as matrices, 
    the corresponding subspaces are the row spaces). In other words, subsp[r][w] 
    is a list of the subspaces D in C with dim(D)=r, |supp(D)|=w, and with trivial
    intersection with C2.
    The optional parameter R allows the user to only compute the relative higher weight 
    spectrum for r in R. If none is given, we compute it for 1 <= r <= dim(C)-dim(C2). 
    This is a version of rhigher_spectrum that requires less memory, at the expense 
    of speed in some cases.

    OUTPUT:

    If subspace_list=False, a dictionary 'spec', whose keys are given by r, and
    whose values are dictionaries, whose keys are the possible weights w, and
    whose values are the number of subcodes D of C with dim(D)=r and |supp(E)|=w.
    If subspace_list=False, a list [spec, subsp], where 'spec' and 'subsp'
    are the dictionaries mentioned above. 

    EXAMPLES:
    
    >> G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
    >> G2 = matrix(GF(2), [G[-1]])
    >> C = LinearCode(G)
    >> C2 = LinearCode(G2)
    >> rhigher_spectrum_low_mem(C, C2)
    {0: {0: 1}, 1: {3: 4, 4: 1, 6: 1}, 2: {5: 2, 6: 2}}
    >> rhigher_spectrum_low_mem(C, C2, verbose=True)
    r = 0
    Spectrum: {0: 1}
    r = 1
    Spectrum: {3: 4, 4: 1, 6: 1}                                                                                                     
    r = 2
    Spectrum: {5: 2, 6: 2}                                                                                                     
    {0: {0: 1}, 1: {3: 4, 4: 1, 6: 1}, 2: {5: 2, 6: 2}}
    >> rhigher_spectrum_low_mem(C, C2, R=[1])
    {1: {3: 4, 4: 1, 6: 1}}
    >> rhigher_spectrum_low_mem(C, C2, R=[2], subspace_list=True)
    [{2: {5: 2, 6: 2}},
     {2: {5: [
    [1 0 0 0 1 1]  [1 0 1 0 0 1]
    [0 1 0 1 1 0], [0 1 1 1 0 0]
    ],
       6: [
    [1 0 1 0 0 1]  [1 0 0 0 1 1]
    [0 1 0 1 1 0], [0 1 1 1 0 0]
    ]}}]

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    H = C.parity_check_matrix()
    G2 = C2.systematic_generator_matrix()
    H2 = C2.parity_check_matrix()
    k2 = C2.dimension()
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    spec = {} 
    if subspace_list:
        subsp = {} 
    if R is None:
        R = range(k - k2 + 1) 
    elif not set(R).issubset(set(range(k - k2 + 1))):
        raise Exception('Invalid range of values for R')
        
    for r in R:
        if verbose:
            print('r =', r)
        if r == 0:
            spectemp = {0: 1}
            spec.update({0: spectemp})
            if subspace_list:
                subsptemp = {0: [matrix(K, [0 for i in range(n)])]}
                subsp = {0: subsptemp}
            if verbose:
                print('Spectrum:', spec[r]) 
            continue  
        spectemp = {}
        if subspace_list:
            subsptemp = {}
        for w in range(r, k + 1):
            y = range(w) # We start with support {1,...,w}
            for s in combinations(y[1:], r - 1): # We assume we have a pivot on the first
                # position, and we choose r - 1 more pivots
                comp = [z for z in y if z not in s and z!=y[0]] # Non pivots
                ss = [y[0]] + list(s) # Pivots
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
                        
                        Mtemp = matrix(K, MM).transpose() * G
                        soptemp = len(matrix_supp(Mtemp))
                        interd = r - (H2 * Mtemp.transpose()).rank()
                        if interd == 0:
                            try:
                                spectemp.update({soptemp: spectemp[soptemp] + 1})
                            except:
                                spectemp.update({soptemp: 1})
                            if subspace_list:
                                try:
                                    subsptemp.update({soptemp: subsptemp[soptemp] + [Mtemp]})
                                except:
                                    subsptemp.update({soptemp: [Mtemp]})
            if verbose:
                spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} # Order by key
                print('Current spectrum:', spectemp, 'w =', w, end='\r') 
        spectemp = {k: v for k, v in sorted(spectemp.items(), key=lambda item: item[0])} 
        spec.update({r: spectemp})
        if subspace_list:
            subsptemp = {k: v for k, v in sorted(subsptemp.items(), key=lambda item: item[0])} 
            subsp.update({r: subsptemp})
        if verbose:
            print('Spectrum:',spec[r], ' '*100) 
    if subspace_list:
        return [spec, subsp]
    else:
        return spec