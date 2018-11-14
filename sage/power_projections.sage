load('utils.sage')

##################################################
#  Input:                                        #
#    * generic monic polynomial g of degree n,   #
#    * generic polynomial a of degree < n,       #
#    * linear form ell: K[y] / (g) -> K given by #
# the vector [ell(1) ell(y) .. ell(y^(n-1))]     #
#    * block dimension parameter 2 <= m < n,     #
#  Output:                                       #
#    * vector [ell(1) ell(a) .. ell(a^(n-1))]    #
##################################################

def power_projections(g, a, ell, m, verbose=False):
    (PolyRing,y) = g.parent().objgen()
    Field = PolyRing.base_ring()
    PolyRingX.<x> = Field[]
    n = g.degree()
    d = ceil(n/m)

    if verbose:
        import time
        t_total = 0.0
        print_separator()
        print " COMPUTING POWER PROJECTIONS"
        print_separator()
        print "Context:"
        print "-- working over", Field
        print "-- degree n =", n
        print "-- block dimension m =", m
        print "-- block degree d =", d
        print "-- input modulus g =", g(y) if (n<30) else "{degree n polynomial}"
        print "-- input polynomial a =", a(y) if (n<30) else "{degree <n polynomial}"
        print "-- input linear form ell =", ell if (n<30) else "{length n vector over the field}"
        print_separator()
        print "Goal: compute power projections [ell(1) ell(a) .. ell(a^(n-1))]"
        print_separator()

    if verbose:
        print "First step:"
        print "  Compute a balanced basis for a and g"
        t_start = time.time()

    B, pow_a = balanced_basis(g, a, m, store_powers=True, verbose=False)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Second step:"
        print "  retrieve the determinant of the balanced basis"
        print "  -> this is a polynomial f(x) of degree n"
        t_start = time.time()

    f = B.determinant()

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Third step:"
        print "  compute ell*P, where P is the numerator associated to B,"
        print "  that is, P(x) = (xId - A)^(-1) Matrix.block([[B], [0]]), where"
        print "  A is the multiplication matrix of a in K[y] / <g> in the"
        print "  monomial basis, and the zero block has dimension (n-m) x m."
        t_start = time.time()

    if verbose:
        print "--Third step, first substep:"
        print "    Compute the row vectors ell * A^k for k in range(d)"
        print "    Algorithm: see Lemma 3.1"
        t_start1 = time.time()

    # TODO prototype: for the moment, naive algorithm computing the matrix A
    A = Matrix(Field,n,n) # initialize S[k] with mxm zero matrix
    # for each 0<=j<n, compute the n coefficients of a*y^j mod g
    for j in range(n):
        coeffs = list(a*y^j % g)
        # pad with zeroes if necessary
        A[:,j] = Matrix(n, 1, coeffs + [0]*(n-len(coeffs)))

    ellAk = copy(ell)
    ellA = [ellAk] # initialize with term for k=0
    for k in range(1,d):
        ellAk = ellAk * A  # now ellAk = ell * A**k
        ellAk.append(ellAk)
    
    

    if verbose:
        t_end1 = time.time()
        print "Computation time: ", t_end1-t_start1
        print "--Third step, second substep:"
        print "    Compute the row vectors ell * A^k for k in range(d)"
        print "    Algorithm: see Lemma 3.1"
        t_start1 = time.time()


    if verbose:
        t_end = time.time()
        print "Computation time (third step): ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    return b
