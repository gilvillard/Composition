load('utils.sage')

##################################################
#  Input:                                        #
#    * generic monic polynomial g of degree n,   #
#    * generic polynomial a of degree < n,       #
#    * nonzero polynomial b of degree < n,       #
#    * block dimension parameter m < n           #
#  Output:                                       #
#    * polynomial h of degree < n such that      #
#       h(a) = b mod g                           #
##################################################

def inverse_composition(g, a, b, m, verbose=False):
    (PolyRing,y) = g.parent().objgen()
    Field = PolyRing.base_ring()
    PolyRingX.<x> = Field[]
    n = g.degree()
    d = ceil(n/m)

    if verbose:
        import time
        t_total = 0.0
        print_separator()
        print " PERFORMING INVERSE COMPOSITION"
        print_separator()
        print "Context:"
        print "-- working over", Field
        print "-- degree n =", n
        print "-- block dimension m =", m
        print "-- block degree d =", d
        print "-- input modulus g =", g(y) if (n<30) else "{degree n polynomial}"
        print "-- input polynomial a =", a(y) if (n<30) else "{degree <n polynomial}"
        print "-- input polynomial b =", b(y) if (n<30) else "{degree <n polynomial}"
        print_separator()
        print "Goal: compute composition inverse h such that h(a) = b mod g"
        print_separator()

    if verbose:
        print "First step:"
        print "  Compute the linearly recurrent sequence S[k] for k in 0 ... 2*d-1 =", 2*d-1, ","
        print "  where S[k] = [I_m  0] A^k [I_m  0].T, with A the multiplication matrix of a"
        print "  in the quotient algebra modulo g, in the monomial basis"
        print "  For a given k in {0...2*d}, this can be seen as the m x m matrix whose"
        print "  row j contains the coefficients of degree 0, 1, ..., m-1 of the "
        print "  polynomial a^k y^j mod g"
        print "Complexity: O(d M(n) + m n)"
        print "  Here, d n log_2(n) loglog_2(n) + m n is about", ceil(RR(d*n*log(n,2)*log(log(n,2),2) + m*n))
        t_start = time.time()

    # polynomial: -g truncated modulo y^m
    g_low = (-g).truncate(m)
    # list of last m coefficients of -g (except leading coeff 1)
    g_high = [-g[n-m+i] for i in range(m)]

    # S = the sequence of matrices, initially only S[0] = identity
    S = [Matrix.identity(Field, m)]
    # f = a^k mod g, initially a
    f = a
    # store the powers of a, used in second step
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
        f = (f * a) % g  # a^k mod g
        pow_a.append(f)

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    ## Naive version
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
        print "  Compute the sequence Sb[k] for k in 0 ... 2*d-1 =", 2*d-1, ","
        print "  where Sb[k] = [I_m  0] A^k vb, with A the multiplication matrix of a"
        print "  in the quotient algebra modulo g as above,"
        print "  and vb the coefficient vector of the right-hand side b(y)."
        print "  Explicitly, Sb[k] is simply formed by the first m"
        print "  coefficients of the product a^k b mod g."
        print "Complexity: O(d M(n))"
        print "  Here, d n log_2(n) loglog_2(n) is about", ceil(RR(d*n*log(n,2)*log(log(n,2),2)))
        t_start = time.time()

    # initially Sb contains the vector for k=0, which is the beginning of the
    # coefficient vector of b, padded with zeroes if necessary
    f = b  # stores the products a^k*b
    Sb = [Matrix(Field, n, 1, list(f) + [0]*(n-f.degree()-1))[:m,0]]
    for k in range(1,2*d):
        f = a*f % g
        Sb.append(Matrix(Field, n, 1, list(f) + [0]*(n-f.degree()-1))[:m,0])

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Third step:"
        # TODO text to be updated
        print "  Compute a column-reduced right matrix generator for the"
        print "  sequence S[k], k in 0 ... 2*d-1. Our assumption (generic input)"
        print "  ensures that this is an m x m matrix of degree d."
        print "  This matrix is a balanced basis as described above."
        print "Complexity: O(MM(m,d) log(d)),"
        print "  where MM(m,d) is the time for multiplication of two"
        print "  univariate polynomial m x m matrices of degree at most d."
        print "  Here, m^omega d log_2(d) is about:"
        print "     ", ceil(RR(m^3*d*log(d,2))), "if omega = 3;"
        print "     ", ceil(RR(m^2.38*d*log(d,2))), "if omega = 2.38;"
        print "     ", ceil(RR(m^2*d*log(d,2))), "if omega = 2."
        t_start = time.time()

    # matrix sequence which will be reconstructed as a fraction of polynomial
    # matrices to give a balanced basis
    # F = sum_{0 <= k < 2d+1} S[k] x^(2d-k)
    S.reverse()
    F = Matrix(PolyRingX, m, m)
    for i in range(m):
        for j in range(m):
            F[i,j] = PolyRingX([S[k][i,j] for k in range(2*d)])

    # matrix sequence which will be reconstructed as a bivariate polynomial
    # of bounded degrees, which solves the inverse composition problem
    Sb.reverse()
    Fb = Matrix(PolyRingX, m, 1)
    for i in range(m):
        Fb[i,0] = PolyRingX([Sb[k][i,0] for k in range(2*d)])

    # perform reconstruction via approximant basis computation, with
    # input matrix [Fb | F | -Identity]
    pmat = Matrix.block([[-Fb, F, -1]])
    appbas = approximant_basis(pmat, 2*d, [2*d] + [0]*(2*m))

    # since the approximant basis 'appbas' is in ordered weak Popov form, we
    # know where to look for the balanced basis and the \tilde{h}
    if appbas[0,0].degree()>0:
        print "Error in inverse composition: there is no bivariate"
        print "polynomial hh(x,y) of the specified degree"
        print "which realizes inverse composition"
    hh = appbas[1:m+1,0]
    B = appbas[1:m+1,1:m+1]

    ## SANITY CHECK
    #val = sum([(hh[i,0](a) * y**i) % g for i in range(m)])
    #if val != b:
    #    print "bivpoly error"

    #for j in range(m):
    #    val = sum([(B[i,j](a) * y**i) % g for i in range(m)])
    #    if val != 0:
    #        print "bb error"

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Fourth step: TODO text to be updated"
        print "  Use balanced basis B and bivariate polynomial hh(x,y)"
        print "  to deduce univariate polynomial h(x)"
        print "Complexity: O~(m^omega d),"
        print "  Here, m^omega d is about:"
        print "     ", ceil(RR(m^3*d)), "if omega = 3;"
        print "     ", ceil(RR(m^2.38*d)), "if omega = 2.38;"
        print "     ", ceil(RR(m^2*d)), "if omega = 2."
        t_start = time.time()

    # use kernel basis (TODO clean this part)
    sysmat = Matrix.block([[hh[1:,0], B[1:,:]]])
    kerbas = approximant_basis(sysmat, m*d, [m*d]+[0]*m)
    coeffs = kerbas[1:,0]
    h = hh[0,0] + B[0,:] * coeffs

    # existence check
    if kerbas[0,0].degree()>0:
        print "Error in inverse composition: no solution"
        print "(linear system fail)"

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    return h, B
