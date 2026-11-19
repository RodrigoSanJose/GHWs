#################################################################################
# Author:
# Rodrigo San-José. Contact: rsanjose@vt.edu
# GitHub repository: https://github.com/RodrigoSanJose/GHWs

# This module provides functions to check whether proposed lower or upper
# bounds for generalized Hamming weights and relative generalized Hamming
# weights are valid.

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
    bch_bound, binary_rank, binary_rows, binary_supp, colwt, information, is_cyclic, matrix_supp,
    standard, subspaces
)


def GHW_bound(C, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether bound is a valid lower or upper bound for the rth GHW of C.
    The argument bound_type must be either 'lower' or 'upper'. If
    bound_type='lower', the algorithm returns False as soon as a subcode with
    cardinality of support lower than bound is found, and returns True as soon
    as the lower bound for the subspaces that have not yet been enumerated
    reaches bound. If bound_type='upper', the algorithm returns True as soon as
    a subcode with cardinality of support at most bound is found, and returns
    False as soon as the lower bound for the subspaces that have not yet been
    enumerated is greater than bound. The optional arguments L and verbose
    follow the same conventions as in GHW.

    OUTPUT:

    True if bound is a valid bound of the specified type, and False otherwise.

    EXAMPLES::

        sage: C = codes.BinaryReedMullerCode(1, 5)
        sage: GHW_bound(C, 2, 24)
        True
        sage: GHW_bound(C, 2, 23)
        True
        sage: GHW_bound(C, 2, 25)
        False
        sage: GHW_bound(C, 2, 24, bound_type='upper')
        True
        sage: GHW_bound(C, 2, 25, bound_type='upper')
        True
        sage: GHW_bound(C, 2, 23, bound_type='upper')
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
    ghwmax = n - k + r
    if bound_type == 'lower':
        if bound <= r:
            return True
        elif bound > ghwmax:
            return False
    else:
        if bound < r:
            return False
        elif bound >= ghwmax:
            return True
    # Only cyclic codes with non-repeated roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    if k > 0 and not gen:
        raise ValueError("L[1] must contain at least one generator matrix")
    if any(j.nrows() != k or j.ncols() != n for j in gen):
        raise Exception('The matrices in L[1] have to be k x n generator matrices of C')
    # Over GF(2), the rows of the matrices in L[1] are used as integers (if they are over GF(2))
    binary = K.order() == 2 and all(j.base_ring() == K for j in gen)
    if binary:
        gen_bits = [binary_rows(j) for j in gen]
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C))
        except:
            ghwlb = ghwlb
    ghwub = ghwmax
    if bound_type == 'lower':
        if ghwlb >= bound:
            return True
    else:
        if ghwlb > bound:
            return False
    if ghwlb >= ghwub:
        if bound_type == 'lower':
            return ghwub >= bound
        return ghwub <= bound
    # For a lower bound we need to reach bound. For an upper bound we need to
    # exceed bound in order to prove that no suitable subcode has support at
    # most bound.
    target = bound if bound_type == 'lower' else bound + 1
    w = r
    while w <= k and ghwlb < ghwub:
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < target:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < target:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= target:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to reach target at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)

        # These are the only matrices that contribute
        gen_reduced = [gen[j] for j in range(len(gen)) if j not in rm]
        if binary: # Over GF(2), their rows as integers (see binary_supp)
            gen_reduced = [gen_bits[j] for j in range(len(gen)) if j not in rm]
        rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
        if binary: # Over GF(2), the positions of the ones in each row (see binary_supp)
            rrefs = [[mat.nonzero_positions_in_row(i) for i in range(r)] for mat in rrefs]
        for y in combinations(range(k), w): # All possible supports of weight w
            # The matrix with the columns of mat in the positions of y (and zero elsewhere),
            # times j, is equal to mat times the rows of j in the positions of y
            if binary:
                gen_rows = gen_reduced
            else:
                gen_rows = [j.matrix_from_rows(y) for j in gen_reduced]
            for mat in rrefs:
                for j in gen_rows:
                    if binary:
                        supptemp = binary_supp(mat, j, y)
                    else:
                        supptemp = len(matrix_supp(mat * j))
                    if supptemp < ghwub:
                        ghwub = supptemp
                        if verbose:
                            print('Subspace with cardinality of support', supptemp, 'found')
                        if bound_type == 'lower' and ghwub < bound:
                            return False
                        elif bound_type == 'upper' and ghwub <= bound:
                            return True
        # Lower bound for the subspaces that have not yet been enumerated
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + max((w + 1) - red[j], 0)
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb >= bound:
            return True
        elif bound_type == 'upper' and ghwlb > bound:
            return False
        w = w + 1
    if bound_type == 'lower':
        return ghwub >= bound
    return ghwub <= bound


def RGHW_bound(C, C2, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether bound is a valid lower or upper bound for the rth RGHW of C
    with respect to C2. The argument bound_type must be either 'lower' or
    'upper'. If bound_type='lower', the algorithm returns False as soon as a
    subcode with cardinality of support lower than bound and trivial
    intersection with C2 is found, and returns True as soon as the lower bound
    for the subspaces that have not yet been enumerated reaches bound. If
    bound_type='upper', the algorithm returns True as soon as a subcode with
    cardinality of support at most bound and trivial intersection with C2 is
    found, and returns False as soon as the lower bound for the subspaces that
    have not yet been enumerated is greater than bound. The optional arguments
    L and verbose follow the same conventions as in RGHW.

    OUTPUT:

    True if bound is a valid bound of the specified type, and False otherwise.

    EXAMPLES::

        sage: G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
        sage: G2 = matrix(GF(2), [G[-1]])
        sage: C = LinearCode(G)
        sage: C2 = LinearCode(G2)
        sage: RGHW_bound(C, C2, 1, 3)
        True
        sage: RGHW_bound(C, C2, 1, 2)
        True
        sage: RGHW_bound(C, C2, 1, 4)
        False
        sage: RGHW_bound(C, C2, 1, 3, bound_type='upper')
        True
        sage: RGHW_bound(C, C2, 1, 4, bound_type='upper')
        True
        sage: RGHW_bound(C, C2, 1, 2, bound_type='upper')
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
    if not K.has_coerce_map_from(C2.base_field()): # C2 can be defined over a subfield of K
        raise Exception('C2 has to be defined over the field of C (or a subfield of it)')
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth RGHW is always between r and the generalized Singleton bound
    ghwmax = n - k + r
    if bound_type == 'lower':
        if bound <= r:
            return True
        elif bound > ghwmax:
            return False
    else:
        if bound < r:
            return False
        elif bound >= ghwmax:
            return True
    # Only cyclic codes with non-repeated roots are considered
    cyc = is_cyclic(C) and is_cyclic(C2) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    if k > 0 and not gen:
        raise ValueError("L[1] must contain at least one generator matrix")
    if any(j.nrows() != k or j.ncols() != n for j in gen):
        raise Exception('The matrices in L[1] have to be k x n generator matrices of C')
    # Over GF(2), the rows of the matrices in L[1] are used as integers (if they are over GF(2))
    binary = K.order() == 2 and all(j.base_ring() == K for j in gen)
    if binary: # The rows of each j in gen and of j * H2^T (for the intersection with C2)
        gen_bits = [[binary_rows(j), binary_rows(j * H2.transpose())] for j in gen]
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C))
        except:
            ghwlb = ghwlb
    ghwub = ghwmax
    if bound_type == 'lower':
        if ghwlb >= bound:
            return True
    else:
        if ghwlb > bound:
            return False
    if ghwlb >= ghwub:
        if bound_type == 'lower':
            return ghwub >= bound
        return ghwub <= bound
    # For a lower bound we need to reach bound. For an upper bound we need to
    # exceed bound in order to prove that no suitable subcode has support at
    # most bound.
    target = bound if bound_type == 'lower' else bound + 1
    w = r
    while w <= k and ghwlb < ghwub:
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < target:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < target:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= target:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to reach target at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)

        # These are the only matrices that contribute
        gen_reduced = [gen[j] for j in range(len(gen)) if j not in rm]
        if binary: # Over GF(2), their rows (and those of j * H2^T) as integers
            gen_reduced = [gen_bits[j] for j in range(len(gen)) if j not in rm]
        rrefs = subspaces(r, w, w, K) # All reduced row echelon forms in the first w columns
        if binary: # Over GF(2), the positions of the ones in each row (see binary_supp)
            rrefs = [[mat.nonzero_positions_in_row(i) for i in range(r)] for mat in rrefs]
        for y in combinations(range(k), w): # All possible supports of weight w
            # The matrix with the columns of mat in the positions of y (and zero elsewhere),
            # times j, is equal to mat times the rows of j in the positions of y
            if binary:
                gen_rows = gen_reduced
            else:
                gen_rows = [j.matrix_from_rows(y) for j in gen_reduced]
            for mat in rrefs:
                for j in gen_rows:
                    if binary: # j contains the rows of j and of j * H2^T as integers
                        supptemp = binary_supp(mat, j[0], y)
                    else:
                        Mtempj = mat * j
                        supptemp = len(matrix_supp(Mtempj))
                    if supptemp < ghwub:
                        if binary:
                            interd = r - binary_rank(mat, j[1], y)
                        else:
                            interd = r - (H2 * Mtempj.transpose()).rank()
                        if interd == 0:
                            if verbose:
                                print('Subspace with cardinality of support', supptemp, 'found')
                            ghwub = supptemp
                            if bound_type == 'lower' and ghwub < bound:
                                return False
                            elif bound_type == 'upper' and ghwub <= bound:
                                return True
        # Lower bound for the subspaces that have not yet been enumerated
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + max((w + 1) - red[j], 0)
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb >= bound:
            return True
        elif bound_type == 'upper' and ghwlb > bound:
            return False
        w = w + 1
    if bound_type == 'lower':
        return ghwub >= bound
    return ghwub <= bound


def GHW_bound_low_mem(C, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether bound is a valid lower or upper bound for the rth GHW of C.
    The argument bound_type must be either 'lower' or 'upper'. If
    bound_type='lower', the algorithm returns False as soon as a subcode with
    cardinality of support lower than bound is found, and returns True as soon
    as the lower bound for the subspaces that have not yet been enumerated
    reaches bound. If bound_type='upper', the algorithm returns True as soon as
    a subcode with cardinality of support at most bound is found, and returns
    False as soon as the lower bound for the subspaces that have not yet been
    enumerated is greater than bound. The optional arguments L and verbose
    follow the same conventions as in GHW_low_mem. This is a version of
    GHW_bound that requires less memory, at the expense of speed in some cases.

    OUTPUT:

    True if bound is a valid bound of the specified type, and False otherwise.

    EXAMPLES::

        sage: C = codes.BinaryReedMullerCode(1, 5)
        sage: GHW_bound_low_mem(C, 2, 24)
        True
        sage: GHW_bound_low_mem(C, 2, 23)
        True
        sage: GHW_bound_low_mem(C, 2, 25)
        False
        sage: GHW_bound_low_mem(C, 2, 24, bound_type='upper')
        True
        sage: GHW_bound_low_mem(C, 2, 25, bound_type='upper')
        True
        sage: GHW_bound_low_mem(C, 2, 23, bound_type='upper')
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
    ghwmax = n - k + r
    if bound_type == 'lower':
        if bound <= r:
            return True
        elif bound > ghwmax:
            return False
    else:
        if bound < r:
            return False
        elif bound >= ghwmax:
            return True
    # Only cyclic codes with non-repeated roots are considered
    cyc = is_cyclic(C) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    if k > 0 and not gen:
        raise ValueError("L[1] must contain at least one generator matrix")
    if any(j.nrows() != k or j.ncols() != n for j in gen):
        raise Exception('The matrices in L[1] have to be k x n generator matrices of C')
    # Over GF(2), the rows of the matrices in L[1] are used as integers (if they are over GF(2))
    binary = K.order() == 2 and all(j.base_ring() == K for j in gen)
    if binary:
        gen_bits = [binary_rows(j) for j in gen]
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C))
        except:
            ghwlb = ghwlb
    ghwub = ghwmax
    if bound_type == 'lower':
        if ghwlb >= bound:
            return True
    else:
        if ghwlb > bound:
            return False
    if ghwlb >= ghwub:
        if bound_type == 'lower':
            return ghwub >= bound
        return ghwub <= bound
    # For a lower bound we need to reach bound. For an upper bound we need to
    # exceed bound in order to prove that no suitable subcode has support at
    # most bound.
    target = bound if bound_type == 'lower' else bound + 1
    w = r
    while w <= k and ghwlb < ghwub:
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < target:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < target:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= target:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to reach target at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)

        # These are the only matrices that contribute
        gen_reduced = [gen[j] for j in range(len(gen)) if j not in rm]
        if binary: # Over GF(2), their rows as integers (see binary_supp)
            gen_reduced = [gen_bits[j] for j in range(len(gen)) if j not in rm]
        standards = [standard(t, r, K) for t in range(r)] # The vectors of the pivot columns
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
                # The r x w matrix with pivots in ss and non-pivot columns i, which is the
                # same for every sup
                MM = []
                for z in range(w):
                    if z in ss: # Pivots
                        MM.append(standards[ss.index(z)])
                    else: # Non pivots
                        MM.append(i[comp.index(z)])
                mat = matrix(K, MM).transpose()
                if binary: # The positions of the ones in each row (see binary_supp)
                    mat = [mat.nonzero_positions_in_row(t) for t in range(r)]
                for sup in combinations(range(k), w): # All possible supports of weight w
                    # The matrix with the columns of mat in the positions of sup (and zero
                    # elsewhere), times j, is equal to mat times the rows of j in the
                    # positions of sup (all of j if w = k, so then j is not copied)
                    for j in gen_reduced:
                        if binary:
                            supptemp = binary_supp(mat, j, sup)
                        else:
                            rows = j if w == k else j.matrix_from_rows(sup)
                            supptemp = len(matrix_supp(mat * rows))
                        if supptemp < ghwub:
                            ghwub = supptemp
                            if verbose:
                                print('Subspace with cardinality of support', supptemp, 'found')
                            if bound_type == 'lower' and ghwub < bound:
                                return False
                            elif bound_type == 'upper' and ghwub <= bound:
                                return True
        # Lower bound for the subspaces that have not yet been enumerated
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + max((w + 1) - red[j], 0)
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb >= bound:
            return True
        elif bound_type == 'upper' and ghwlb > bound:
            return False
        w = w + 1
    if bound_type == 'lower':
        return ghwub >= bound
    return ghwub <= bound


def RGHW_bound_low_mem(C, C2, r, bound, bound_type='lower', L=None, verbose=False):
    r"""
    Checks whether bound is a valid lower or upper bound for the rth RGHW of C
    with respect to C2. The argument bound_type must be either 'lower' or
    'upper'. If bound_type='lower', the algorithm returns False as soon as a
    subcode with cardinality of support lower than bound and trivial
    intersection with C2 is found, and returns True as soon as the lower bound
    for the subspaces that have not yet been enumerated reaches bound. If
    bound_type='upper', the algorithm returns True as soon as a subcode with
    cardinality of support at most bound and trivial intersection with C2 is
    found, and returns False as soon as the lower bound for the subspaces that
    have not yet been enumerated is greater than bound. The optional arguments
    L and verbose follow the same conventions as in RGHW_low_mem. This is a
    version of RGHW_bound that requires less memory, at the expense of speed in
    some cases.

    OUTPUT:

    True if bound is a valid bound of the specified type, and False otherwise.

    EXAMPLES::

        sage: G = matrix(GF(2), [(1, 0, 0, 0, 1, 1), (0, 1, 0, 1, 1, 0), (0, 0, 1, 0, 1, 0)])
        sage: G2 = matrix(GF(2), [G[-1]])
        sage: C = LinearCode(G)
        sage: C2 = LinearCode(G2)
        sage: RGHW_bound_low_mem(C, C2, 1, 3)
        True
        sage: RGHW_bound_low_mem(C, C2, 1, 2)
        True
        sage: RGHW_bound_low_mem(C, C2, 1, 4)
        False
        sage: RGHW_bound_low_mem(C, C2, 1, 3, bound_type='upper')
        True
        sage: RGHW_bound_low_mem(C, C2, 1, 4, bound_type='upper')
        True
        sage: RGHW_bound_low_mem(C, C2, 1, 2, bound_type='upper')
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
    if not K.has_coerce_map_from(C2.base_field()): # C2 can be defined over a subfield of K
        raise Exception('C2 has to be defined over the field of C (or a subfield of it)')
    if not H * G2.transpose() == 0:
        raise Exception('C2 is not contained in C')
    elif C.dimension() == C2.dimension():
        raise Exception('C cannot be equal to C2')
    if bound_type not in ['lower', 'upper']:
        raise Exception("bound_type has to be either 'lower' or 'upper'")
    # The rth RGHW is always between r and the generalized Singleton bound
    ghwmax = n - k + r
    if bound_type == 'lower':
        if bound <= r:
            return True
        elif bound > ghwmax:
            return False
    else:
        if bound < r:
            return False
        elif bound >= ghwmax:
            return True
    # Only cyclic codes with non-repeated roots are considered
    cyc = is_cyclic(C) and is_cyclic(C2) and list(G.pivots()) == srange(k)
    if L is None:
        if cyc:
            L = [[i + 1 for i in range(k)], [G], [0]]
        else:
            L = information(G)
    [inf, gen, red] = L
    if k > 0 and not gen:
        raise ValueError("L[1] must contain at least one generator matrix")
    if any(j.nrows() != k or j.ncols() != n for j in gen):
        raise Exception('The matrices in L[1] have to be k x n generator matrices of C')
    # Over GF(2), the rows of the matrices in L[1] are used as integers (if they are over GF(2))
    binary = K.order() == 2 and all(j.base_ring() == K for j in gen)
    if binary: # The rows of each j in gen and of j * H2^T (for the intersection with C2)
        gen_bits = [[binary_rows(j), binary_rows(j * H2.transpose())] for j in gen]
    ghwlb = r
    if cyc:
        try:
            ghwlb = max(ghwlb, bch_bound(C))
        except:
            ghwlb = ghwlb
    ghwub = ghwmax
    if bound_type == 'lower':
        if ghwlb >= bound:
            return True
    else:
        if ghwlb > bound:
            return False
    if ghwlb >= ghwub:
        if bound_type == 'lower':
            return ghwub >= bound
        return ghwub <= bound
    # For a lower bound we need to reach bound. For an upper bound we need to
    # exceed bound in order to prove that no suitable subcode has support at
    # most bound.
    target = bound if bound_type == 'lower' else bound + 1
    w = r
    while w <= k and ghwlb < ghwub:
        # Computation of w0, the expected w to finish
        rm = []
        if cyc:
            w0 = w
            while ceil((w0 + 1) * n / k) < target:
                w0 = w0 + 1
        else:
            w0 = w - 1
            ghwlbtemp = 0
            while ghwlbtemp < target:
                ghwlbtemp = 0
                w0 = w0 + 1
                for j in range(len(gen)):
                    ghwlbtemp = ghwlbtemp + max((w0 + 1) - red[j], 0)
                    if ghwlbtemp >= target:
                        if w0 == w:
                            # We store in rm the indices of the matrices that are not necessary
                            # to reach target at the end of this iteration (if any)
                            rm = rm + srange(j + 1, len(gen))
                        break
        if verbose:
            print('Lower:', ghwlb, 'Upper:', ghwub, 'Support:', w, 'Expected:', w0)

        # These are the only matrices that contribute
        gen_reduced = [gen[j] for j in range(len(gen)) if j not in rm]
        if binary: # Over GF(2), their rows (and those of j * H2^T) as integers
            gen_reduced = [gen_bits[j] for j in range(len(gen)) if j not in rm]
        standards = [standard(t, r, K) for t in range(r)] # The vectors of the pivot columns
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
                # The r x w matrix with pivots in ss and non-pivot columns i, which is the
                # same for every sup
                MM = []
                for z in range(w):
                    if z in ss: # Pivots
                        MM.append(standards[ss.index(z)])
                    else: # Non pivots
                        MM.append(i[comp.index(z)])
                mat = matrix(K, MM).transpose()
                if binary: # The positions of the ones in each row (see binary_supp)
                    mat = [mat.nonzero_positions_in_row(t) for t in range(r)]
                for sup in combinations(range(k), w): # All possible supports of weight w
                    # The matrix with the columns of mat in the positions of sup (and zero
                    # elsewhere), times j, is equal to mat times the rows of j in the
                    # positions of sup (all of j if w = k, so then j is not copied)
                    for j in gen_reduced:
                        if binary: # j contains the rows of j and of j * H2^T as integers
                            supptemp = binary_supp(mat, j[0], sup)
                        else:
                            rows = j if w == k else j.matrix_from_rows(sup)
                            Mtempj = mat * rows
                            supptemp = len(matrix_supp(Mtempj))
                        if supptemp < ghwub:
                            if binary:
                                interd = r - binary_rank(mat, j[1], sup)
                            else:
                                interd = r - (H2 * Mtempj.transpose()).rank()
                            if interd == 0:
                                if verbose:
                                    print('Subspace with cardinality of support', supptemp, 'found')
                                ghwub = supptemp
                                if bound_type == 'lower' and ghwub < bound:
                                    return False
                                elif bound_type == 'upper' and ghwub <= bound:
                                    return True
        # Lower bound for the subspaces that have not yet been enumerated
        ghwlbtemp = 0
        for j in range(len(gen_reduced)):
            if cyc:
                ghwlbtemp = ceil((w + 1) * n / k)
            else:
                ghwlbtemp = ghwlbtemp + max((w + 1) - red[j], 0)
        ghwlb = max(ghwlb, ghwlbtemp)
        if bound_type == 'lower' and ghwlb >= bound:
            return True
        elif bound_type == 'upper' and ghwlb > bound:
            return False
        w = w + 1
    if bound_type == 'lower':
        return ghwub >= bound
    return ghwub <= bound
