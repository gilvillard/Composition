############################
#  SOLVING LINEAR SYSTEMS  #
############################

# Below, in cost bounds,
# Newton(m,d) = MM(m,d) + MM(m,d/2) + MM(m,d/4) + ... + MM(m,1),
# where MM is the cost for polynomial matrix multiplication
# in dimension m and degree d

# Truncated inverse X = A^{-1} mod x^d
def inverse_truncated(A, d):
    ##Input:
    #    * m x m polynomial matrix A, with A(0) invertible
    #    * positive integer d (order of truncation)
    ##Output:
    #    * the m x m polynomial matrix X such that X = A^(-1) mod d
    ##Algorithm: Newton iteration
    ##Complexity: O(Newton(m,d))
    Id = Matrix.identity(A.base_ring(), A.nrows())
    X = A(0).inverse() # B = A.inverse() mod x^1
    for i in range(1,1+d.nbits()):
        #Newton iteration: B = (2I - BA)B mod x^(2^i)
        X = (2*Id - X*A) * X
        matrix_truncate(X, 2**i)
    # we might have gone a little too far in the truncation order
    # --> truncate to respect the algorithm specification
    matrix_truncate(X, d)
    return X

def system_solve_expansion(A, B, d):
    ##Input:
    #    * m x m polynomial matrix A, with A(0) invertible
    #    * m x n polynomial matrix B
    #    * positive integer d (order of truncation)
    ##Output:
    #    * the m x n polynomial matrix X such that A*X = B mod x^d
    ##Algorithm: doing d / deg(A) steps at precision deg(A)
    ##Complexity (assuming deg(B) in O(m/n deg(A)):
    #       O( Newton(m,deg(A)) + d/deg(A) MM(m, deg(A)) ) 
    m = A.ncols()
    dA = 1+A.degree()
    Ai = inverse_truncated(A, dA)
    BB = copy(B)
    X = Matrix(A.base_ring(), m, B.ncols())
    for k in range(0,(d/dA).ceil()):
        # compute XX = Ai * BB mod x^dA
        XX = Ai * BB   # Note: asymptotically fast computation would require partial linearization on BB
        matrix_truncate(XX, dA)
        # compute BB = (BB - A*XX) x^(-dA)
        BB = BB - A*XX
        matrix_shift(BB, -dA)
        # update X = X + x^(k*dA) * XX
        matrix_shift(XX, k*dA)
        X = X + XX
    matrix_truncate(X, d)
    return X

def system_solve(A, B):
    ##Input:
    #    * m x m polynomial matrix A, with A(0) invertible
    #    * m x 1 polynomial vector B
    ##Output:
    #    * the couple (X,f) where X is an m x 1 polynomial vector
    #      and f is the minimal monic univariate polynomial 
    #      such that A*X = f*B
    ##Algorithm:
    #   compute an expansion of f^(-1) X at sufficiently large order,
    #   then reconstruct the denominator f and the numerator X.
    ##Complexity (assuming deg(B) <=  m deg(A)):
    #       O( Newton(m,deg(A)) + m MM(m, deg(A)) ) 
    var = A.base_ring().gen()   # the variable
    m = A.ncols()
    dA = m * A.degree()  # bound on deg(f)
    dB = B.degree()  # dA+dB is a bound on deg(X)
    # find X = A^{-1} B mod x^d for sufficiently large d
    d = 2*dA+dB+1
    X = system_solve_expansion(A, B, d)
    f = 1
    # find the common denominator f
    for k in range(m):
        # reconstruct fraction X[k,0] = N/D mod x^d,
        # with deg(D) <= dA, deg(N) <= dA+dB
        N,D = X[k,0].rational_reconstruct(var^d, dA+dB, dA)
        f = lcm(f,D)
        if f.degree() == dA: # we have found f, early termination
            break
    # build corresponding numerators
    for k in range(m):
        X[k,0] = (f * X[k,0]).truncate(dA+dB+1)
    return (X,f)
