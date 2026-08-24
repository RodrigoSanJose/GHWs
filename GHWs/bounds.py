#################################################################################
# Author:
# Rodrigo San-José. Contact: rsanjose@vt.edu
# GitHub repository: https://github.com/RodrigoSanJose/GHWs

# This module provides functions to check whether known lower or upper bounds
# for generalized Hamming weights and relative generalized Hamming weights are
# attained.

# Required imports
# Standard library
from itertools import combinations

# Sage imports
from sage.arith.srange import srange
from sage.functions.other import ceil
from sage.categories.sets_cat import cartesian_product
from sage.matrix.constructor import matrix

# GHWs imports
from .core import (
    bch_bound, colwt, information, is_cyclic, matrix_supp, standard, subspaces
)


def GHW_bound(C, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether a known lower or upper bound for the rth GHW of C is
    attained. The argument bound_type must be either 'lower' or 'upper'. The
    bound is assumed to be valid. If bound_type='lower', the algorithm returns
    True as soon as a subcode with cardinality of support equal to bound is
    found. If bound_type='upper', the algorithm returns True as soon as the
    lower bound obtained by the algorithm reaches bound. The optional arguments
    L and verbose follow the same conventions as in GHW.

    OUTPUT:

    True if the bound is attained, and False otherwise.

    EXAMPLES::

        sage: C = codes.BinaryReedMullerCode(1, 5)
        sage: GHW_bound(C, 2, 24)
        True
        sage: GHW_bound(C, 2, 23)
        False
        sage: GHW_bound(C, 2, 24, bound_type='upper')
        True
        sage: GHW_bound(C, 2, 25, bound_type='upper')
        False

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception('Invalid value of r')
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth GHW is always between r and the generalized Singleton bound
    if bound not in range(r, n - k + r + 1):
        return False
    # Only cyclic codes with non-repeated roots are considered
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
    if bound_type == 'lower':
        ghwlb = max(ghwlb, bound)
        if ghwlb > bound:
            return False
    else:
        ghwub = min(ghwub, bound)
        if ghwlb > bound:
            return False
    if ghwlb == ghwub:
        return ghwub == bound
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
                    supptemp = len(matrix_supp(Mtemp * j))
                    if supptemp < ghwub:
                        ghwub = supptemp
                        if verbose:
                            print('Subspace with cardinality of support', supptemp, 'found')
                        if bound_type == 'lower':
                            if ghwub == bound:
                                return True
                            elif ghwub < bound:
                                return False
                        else:
                            return False
        # Lower bound calculations
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb > bound:
            return False
        elif bound_type == 'upper' and ghwlb >= bound:
            return ghwlb == bound
        w = w + 1
    return ghwlb == ghwub and ghwub == bound


def RGHW_bound(C, C2, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether a known lower or upper bound for the rth RGHW of C with
    respect to C2 is attained. The argument bound_type must be either 'lower' or
    'upper'. The bound is assumed to be valid. If bound_type='lower', the
    algorithm returns True as soon as a subcode with cardinality of support
    equal to bound and trivial intersection with C2 is found. If
    bound_type='upper', the algorithm returns True as soon as the lower bound
    obtained by the algorithm reaches bound. The optional arguments L and
    verbose follow the same conventions as in RGHW.

    OUTPUT:

    True if the bound is attained, and False otherwise.

    EXAMPLES::

        sage: G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
        sage: G2 = matrix(GF(2), [G[-1]])
        sage: C = LinearCode(G)
        sage: C2 = LinearCode(G2)
        sage: RGHW_bound(C, C2, 2, 5)
        True
        sage: RGHW_bound(C, C2, 2, 4)
        False
        sage: RGHW_bound(C, C2, 2, 5, bound_type='upper')
        True
        sage: RGHW_bound(C, C2, 2, 6, bound_type='upper')
        False

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
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth RGHW is always between r and the generalized Singleton bound
    if bound not in range(r, n - k + r + 1):
        return False
    # Only cyclic codes with non-repeated roots are considered
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
    if bound_type == 'lower':
        ghwlb = max(ghwlb, bound)
        if ghwlb > bound:
            return False
    else:
        ghwub = min(ghwub, bound)
        if ghwlb > bound:
            return False
    if ghwlb == ghwub:
        return ghwub == bound
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
                    supptemp = len(matrix_supp(Mtempj))
                    if supptemp < ghwub:
                        interd = r - (H2 * Mtempj.transpose()).rank()
                        if interd == 0:
                            if verbose:
                                print('Subspace with cardinality of support', supptemp, 'found')
                            ghwub = supptemp
                            if bound_type == 'lower':
                                if ghwub == bound:
                                    return True
                                elif ghwub < bound:
                                    return False
                            else:
                                return False
        # Lower bound calculations
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb > bound:
            return False
        elif bound_type == 'upper' and ghwlb >= bound:
            return ghwlb == bound
        w = w + 1
    return ghwlb == ghwub and ghwub == bound


def GHW_bound_low_mem(C, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether a known lower or upper bound for the rth GHW of C is
    attained. The argument bound_type must be either 'lower' or 'upper'. The
    bound is assumed to be valid. If bound_type='lower', the algorithm returns
    True as soon as a subcode with cardinality of support equal to bound is
    found. If bound_type='upper', the algorithm returns True as soon as the
    lower bound obtained by the algorithm reaches bound. The optional arguments
    L and verbose follow the same conventions as in GHW_low_mem. This is a
    version of GHW_bound that requires less memory, at the expense of speed in
    some cases.

    OUTPUT:

    True if the bound is attained, and False otherwise.

    EXAMPLES::

        sage: C = codes.BinaryReedMullerCode(1, 5)
        sage: GHW_bound_low_mem(C, 2, 24)
        True
        sage: GHW_bound_low_mem(C, 2, 23)
        False
        sage: GHW_bound_low_mem(C, 2, 24, bound_type='upper')
        True
        sage: GHW_bound_low_mem(C, 2, 25, bound_type='upper')
        False

    """
    K = C.base_field()
    k = C.dimension()
    n = C.length()
    G = C.systematic_generator_matrix()
    if r not in range(1, k + 1):
        raise Exception('Invalid value of r')
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth GHW is always between r and the generalized Singleton bound
    if bound not in range(r, n - k + r + 1):
        return False
    # Only cyclic codes with non-repeated roots are considered
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
    if bound_type == 'lower':
        ghwlb = max(ghwlb, bound)
        if ghwlb > bound:
            return False
    else:
        ghwub = min(ghwub, bound)
        if ghwlb > bound:
            return False
    if ghwlb == ghwub:
        return ghwub == bound
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
                        supptemp = len(matrix_supp(Mtemp * j))
                        if supptemp < ghwub:
                            ghwub = supptemp
                            if verbose:
                                print('Subspace with cardinality of support', supptemp, 'found')
                            if bound_type == 'lower':
                                if ghwub == bound:
                                    return True
                                elif ghwub < bound:
                                    return False
                            else:
                                return False
        # Lower bound calculations
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb > bound:
            return False
        elif bound_type == 'upper' and ghwlb >= bound:
            return ghwlb == bound
        w = w + 1
    return ghwlb == ghwub and ghwub == bound


def RGHW_bound_low_mem(C, C2, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether a known lower or upper bound for the rth RGHW of C with
    respect to C2 is attained. The argument bound_type must be either 'lower' or
    'upper'. The bound is assumed to be valid. If bound_type='lower', the
    algorithm returns True as soon as a subcode with cardinality of support
    equal to bound and trivial intersection with C2 is found. If
    bound_type='upper', the algorithm returns True as soon as the lower bound
    obtained by the algorithm reaches bound. The optional arguments L and
    verbose follow the same conventions as in RGHW_low_mem. This is a version of
    RGHW_bound that requires less memory, at the expense of speed in some cases.

    OUTPUT:

    True if the bound is attained, and False otherwise.

    EXAMPLES::

        sage: G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
        sage: G2 = matrix(GF(2), [G[-1]])
        sage: C = LinearCode(G)
        sage: C2 = LinearCode(G2)
        sage: RGHW_bound_low_mem(C, C2, 2, 5)
        True
        sage: RGHW_bound_low_mem(C, C2, 2, 4)
        False
        sage: RGHW_bound_low_mem(C, C2, 2, 5, bound_type='upper')
        True
        sage: RGHW_bound_low_mem(C, C2, 2, 6, bound_type='upper')
        False

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
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth RGHW is always between r and the generalized Singleton bound
    if bound not in range(r, n - k + r + 1):
        return False
    # Only cyclic codes with non-repeated roots are considered
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
    if bound_type == 'lower':
        ghwlb = max(ghwlb, bound)
        if ghwlb > bound:
            return False
    else:
        ghwub = min(ghwub, bound)
        if ghwlb > bound:
            return False
    if ghwlb == ghwub:
        return ghwub == bound
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
                        supptemp = len(matrix_supp(Mtempj))
                        if supptemp < ghwub:
                            interd = r - (H2 * Mtempj.transpose()).rank()
                            if interd == 0:
                                if verbose:
                                    print('Subspace with cardinality of support', supptemp, 'found')
                                ghwub = supptemp
                                if bound_type == 'lower':
                                    if ghwub == bound:
                                        return True
                                    elif ghwub < bound:
                                        return False
                                else:
                                    return False
        # Lower bound calculations
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + (w + 1) - red[j]
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb > bound:
            return False
        elif bound_type == 'upper' and ghwlb >= bound:
            return ghwlb == bound
        w = w + 1
    return ghwlb == ghwub and ghwub == bound
