load('utils.sage')

##################################################
#  Input:                                        #
#    * generic monic polynomial g of degree n,   #
#    * generic polynomial a of degree < n,       #
#    * block dimension parameter 2 <= m < n,     #
#  Output (default, if store_powers=False):      #
#    * a balanced basis for a and g              #
#  Output (if store_powers=True):                #
#    * a balanced basis for a and g              #
#    * the list of modular powers                #
#        a^k mod g for k in 0 ... 2*ceil(n/m)    #
##################################################

def balanced_basis(g, a, m, store_powers=False, verbose=False):
    (PolyRing,y) = g.parent().objgen()
    Field = PolyRing.base_ring()
    PolyRingX.<x> = Field[]
    n = g.degree()
    d = ceil(n/m)

    if verbose:
        import time
        t_total = 0.0
        print_separator()
        print " COMPUTING BALANCED BASIS"
        print_separator()
        print "Context:"
        print "-- working over", Field
        print "-- degree n =", n
        print "-- block dimension m =", m
        print "-- block degree d =", d
        print "-- input modulus g =", g(y) if (n<30) else "{degree n polynomial}"
        print "-- input polynomial a =", a(y) if (n<30) else "{degree <n polynomial}"
        print_separator()
        print "Goal: find balanced basis of a and g, that is:"
        print "  a", m, "x", m, "univariate polynomial matrix of degree a most", d, ","
        print "  whose columns form a basis of the K[x]-module of rank", m, "formed"
        print "  by the polynomials f in <a(y), x-g(y)> with deg_y(f) <", m
        print_separator()

    if verbose:
        print "First step:"
        print "  Compute the linearly recurrent sequence S[k] for k in 0 ... 2*d-1 =", 2*d-1, ","
        print "  where S[k] = [I_m  0] A^k [I_m  0].T, with A the multiplication matrix of a"
        print "  in the quotient algebra modulo g, in the monomial basis"
        print "  For a given k in {0...2*d}, this can be seen as the m x m matrix whose"
        print "  row j contains the coefficients of degree 0, 1, ..., m-1 of the "
        print "  polynomial a^k y^j mod g"
        t_start = time.time()

    # polynomial: -g truncated modulo y^m
    g_low = (-g).truncate(m)
    # list of last m coefficients of -g (except leading coeff 1)
    g_high = [-g[n-m+i] for i in range(m)]

    # S = the sequence of matrices, initially only S[0] = identity
    S = [Matrix.identity(Field, m)]
    # f = a^k mod g, initially a
    f = a
    if store_powers:
        pow_a = [PolyRing(1), f]

    for k in range(1, 2*d):
        # initialize S[k] with mxm zero matrix
        S_k = Matrix(Field,m,m)
        # for each 0<=j<m, compute first m coefficients of a*y^j mod g

        # polynomial: a^k mod g truncated modulo y^m
        f_low = f.truncate(m)
        # list of last m coefficients of a^k mod g
        f_high = [f[n-m+i] for i in range(m)]

        for j in range(m):
            # at this point, f_low = first m coefficients of y^j a^k mod g
            for i in range(m):
                S_k[i,j] = f_low[i]
            # compute the new f_low: first m coefficients of y^{j+1} a^k mod g
            # this is (y * f_low) % y^m  +  f_high[m-1] * g_low
            f_low = f_low.shift(1).truncate(m); # f_low = (y * f_low) % y^m
            lt = f_high[m-1]; # leading term of f_high
            f_low += lt * g_low;
            # also update f_high  (loop goes from m-1 downto 1)
            for i in range(m-1,0,-1):
                f_high[i] = f_high[i-1] + lt * g_high[i];
        S.append(S_k)
        if k != 2*d-1:
            f = (f * a) % g  # a^k mod g
            if store_powers:
                pow_a.append(f)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    ## Naive version, much slower
    #SS = []  # sequence of matrices
    #ff = 1   # powers a^k mod g
    #for k in range(2*d):
    #    # initialize S[k] with mxm zero matrix
    #    SS_k = Matrix(Field,m,m)
    #    # for each 0<=j<m, compute first m coefficients of a*y^j mod g
    #    for j in range(m):
    #        coeffs = list(ff*y^j % g)
    #        # take first m coeffs, padding with zeroes if necessary
    #        SS_k[:,j] = Matrix(m, 1, coeffs[:m] + [0]*(m-len(coeffs)))
    #    SS.append(SS_k)
    #    ff = (ff * a) % g

    if verbose:
        print "Second step:"
        print "  Compute a column-reduced right matrix generator for the"
        print "  sequence S[k], k in 0 ... 2*d-1. Our assumption (generic input)"
        print "  ensures that this is an m x m matrix of degree d."
        print "  This matrix is a balanced basis as described above."
        t_start = time.time()

    # reconstruct sequence as a fraction of polynomial matrices
    # F = sum_{0 <= k < 2d+1} S[k] x^(2d-k)
    S.reverse()
    F = Matrix(PolyRingX, m, m)
    for i in range(m):
        for j in range(m):
            F[i,j] = PolyRingX([S[k][i,j] for k in range(2*d)])

    # do reconstruction via approximant basis computation, with
    # input matrix [F | -Identity]
    pmat = Matrix.block( [[F, -1]] )
    appbas = approximant_basis(pmat, 2*d)

    # since the approximant basis 'appbas' is in ordered weak Popov form, we
    # know that the fraction is in its leftmost rows: the sought generator is
    # the denominator, which is the leading principal mxm submatrix
    B = appbas[:m,:m]

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if store_powers:
        return B, pow_a
    else:
        return B
