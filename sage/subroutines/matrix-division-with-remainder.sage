###############################################
#  POLYNOMIAL MATRIX DIVISION WITH REMAINDER  #
###############################################

def matrix_quo_rem(A, B):
    ##Input:
    #   * an m x n univariate polynomial matrix A,
    #   * an m x m univariate polynomial matrix B,
    ##Requirement: B is in row reduced form
    ##Output:
    #   * the unique couple (Q,R) of m x n polynomial matrices such that
    #     A = BQ + R, with rdeg(R) < rdeg(A) entrywise
    ##Reference: we follow the folklore algorithm which generalizes the fast
    # division of univariate polynomials; precisely we implement [Neiger-Vu,
    # ISSAC 2017, Algorithm 1], adapted to the case of left-division
    m = B.nrows()
    n = A.ncols()
    PolyRing = A.base_ring()
    # Step 0: find parameter d  (delta in above reference)
    rdegA = A.row_degrees() # zero rows of A --> entries -1 in rdegA
    rdeg = B.row_degrees()  # all non-negative since row reduced
    d = max([rdegA[i]-rdeg[i]+1 for i in range(m)])
    if d<=0: # A already reduced modulo B, quotient is zero
        return (Matrix.zero(PolyRing,m,n), A)
    # Step 1: reverse input matrices
    # Brev = diag(x^(rdeg[i])) B(1/x)
    Brev = Matrix(PolyRing, m, m, [[B[i,j].reverse(rdeg[i]) \
                        for j in range(m)] for i in range(m)])
    # Arev = diag(x^(d+rdeg[i]-1)) A(1/x)
    Arev = Matrix(PolyRing, m, n, [[A[i,j].reverse(d+rdeg[i]-1) \
                        for j in range(n)] for i in range(m)])
    # Step 2: compute quotient
    # Qrev = Brev^{-1} Arev mod x^d 
    Qrev = system_solve_expansion(Brev, Arev, d)
    # quotient is the reverse Q = x^(d-1) Qrev(1/x)
    Q = Matrix(PolyRing, m, n, [[Qrev[i,j].reverse(d-1) \
                    for j in range(n)] for i in range(m)])
    # Step 3: deduce remainder and return
    R = A - B*Q
    return Q,R

############################
#  SOLVING LINEAR SYSTEMS  #
############################

# (used above for matrix division with remainder)

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
