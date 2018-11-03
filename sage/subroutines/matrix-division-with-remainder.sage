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
