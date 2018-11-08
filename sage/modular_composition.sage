load('utils.sage')

##################################################
#  Input:                                        #
#    * generic monic polynomial g of degree n,   #
#    * generic polynomial a of degree < n,       #
#    * generic polynomial h of degree < n,       #
#    * block dimension parameter 2 <= m < n,     #
#  Output:                                       #
#    * h(a) mod g                                #
##################################################

def modular_composition(g, a, h, m, verbose=False):
    (PolyRing,y) = g.parent().objgen()
    Field = PolyRing.base_ring()
    PolyRingX.<x> = Field[]
    n = g.degree()
    d = ceil(n/m)

    if verbose:
        import time
        t_total = 0.0
        print_separator()
        print " PERFORMING MODULAR COMPOSITION"
        print_separator()
        print "Context:"
        print "-- working over", Field
        print "-- degree n =", n
        print "-- block dimension m =", m
        print "-- block degree d =", d
        print "-- input modulus g =", g(y) if (n<30) else "{degree n polynomial}"
        print "-- input polynomial a =", a(y) if (n<30) else "{degree <n polynomial}"
        print "-- input polynomial h =", h(y) if (n<30) else "{degree <n polynomial}"
        print_separator()
        print "Goal: compute modular composition h(a) mod g"
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
        print "  find the remainder of h(x) modulo the balanced basis"
        print "  -> this gives a vector of length m and degree < d"
        print "  seen as a bivariate polynomial of y-degree < m"
        print "  and x-degree < d"
        t_start = time.time()

    # build column vector [h 0 .. 0].T, and reduce it modulo B
    A = Matrix(PolyRingX, m, 1, [h] + [0]*(m-1))
    Q,hh = matrix_quo_rem(A,B)
    # we will not use A and Q anymore

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    if verbose:
        print "Third step:"
        print "  evaluate the bivariate polynomial hh(x,y) obtained at the"
        print "  second step:  h(a) == hh(a,y) mod g"
        t_start = time.time()

    # hh is a vector of length m and degree <d in x
    # rewrite it as a vector of length d and degree <m in y
    H = [sum([hh[i,0][k]*y^i for i in range(m)]) for k in range(d)]

    # compose naively; remember pow_a[k]  ==  a^k % g
    b = sum([H[k] * pow_a[k] for k in range(d)])
    b = b % g

    if verbose:
        t_end = time.time()
        print "Computation time: ", t_end-t_start
        t_total += t_end-t_start
        print_separator()

    return b
