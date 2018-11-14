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

    B = balanced_basis(g, a, m, verbose=False)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Second step:"
        print "  Retrieve the determinant of the balanced basis B(x)"
        print "  -> this is a polynomial f(x) of degree n, resultant"
        print "  in x of x-a(y) and g(y)."
        print "  Also retrieve the quotient vector v_q(x) defined as"
        print "  v_q(x) = B(x)^(-1) [f(x) 0 .. 0].T"
        print "  -> this is a mx1 polynomial vector of degree <= n"
        t_start = time.time()

    f = B.determinant()
    v_f = Matrix(PolyRingX, m, 1, [f]+[0]*(m-1))
    v_q,v_r = matrix_quo_rem(v_f, B)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Third step:"
        print "  compute ell*P, where P is the numerator associated to B,"
        print "  that is, P(x) = C(x) B(x), where C consists of the first m"
        print "  columns of (xId - A)^(-1), where A is the multiplication"
        print "  matrix of a in K[y] / <g> in the monomial basis."
        print "  P is an n x m matrix of degree less than d."
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
        A[:,j] = Matrix(Field, n, 1, coeffs + [0]*(n-len(coeffs)))

    ellAk = Matrix(Field, 1, n, ell)
    ellA = [ellAk] # initialize with term for k=0
    for k in range(1,d+1):
        ellAk = ellAk * A  # now ellAk = ell * A**k
        ellA.append(ellAk)

    if verbose:
        t_end1 = time.time()
        print "--Computation time: ", t_end1-t_start1
        print "--Third step, second substep:"
        print "    deduce ell*P by polynomial vector-matrix multiplication:"
        print "    for v(x) the first m entries of the vector"
        print "    sum_{0 <= k <= d} ell * A^k x^{-k-1}, we have"
        print "    ell*P = terms of degree >= 0 of v(x) * B(x)"
        t_start1 = time.time()

    # compute v(x) = first m entries of the vector
    #     sum_{0 <= k <= d} ell * A^k x^{d-k}"
    # (note the nonnegative exponent d-k, rather than -k-1,
    # used for convenience in the code)
    # --> we will look at terms of degree >= d+1
    ellA.reverse()
    v = sum([ellA[k][0,:m] * x**k for k in range(d+1)])
    ellP = v * B
    matrix_shift(ellP, -d-1)

    if verbose:
        t_end1 = time.time()
        print "--Computation time: ", t_end1-t_start1
        t_end = time.time()
        print "Computation time (total third step): ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Fourth step:"
        print "  Compute the numerator num = ell * P * v_q, and deduce the" 
        print "  power projections by expanding the fraction num(x) / f(x)."
        t_start = time.time()

    num = ellP * v_q

    # define power series ring, at precision n since we want n terms ell(1),...,ell(a^{n-1})
    PowSer.<tt> = PowerSeriesRing(Field, default_prec=n)
    series = num(tt) / f(tt)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    return series
