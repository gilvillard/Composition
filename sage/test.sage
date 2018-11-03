load('utils.sage')

print "#################################"
print "# Testing LINEAR SYSTEM SOLVING #"
print "#################################"

def test_inverse_truncated():
    p = next_prime(1000000, proof=True)  # using big prime to ensure Prob(A(0) singular) ~ 0
    PolyRing.<x> = GF(p)[]
    m = ZZ.random_element(3,10)
    deg = ZZ.random_element(20,40)
    d = ZZ.random_element(5,40)
    A = Matrix.random(PolyRing, m, m, degree=deg)
    X = inverse_truncated(A, d)
    B = A*X  # should be identity mod x^d
    matrix_truncate(B, d)
    return B == Matrix.identity(PolyRing, m)

print "Testing truncated inverse via Newton iteration..."
correct = all([test_inverse_truncated() for i in range(50)])
print " -->", "correct" if correct else "wrong"

def test_system_solve_expansion():
    p = next_prime(1000000, proof=True)  # using big prime to ensure Prob(A(0) singular) ~ 0
    PolyRing.<x> = GF(p)[]
    m = ZZ.random_element(3,10)
    n = ZZ.random_element(1,5)
    deg = ZZ.random_element(20,40)
    d = floor(m/n * deg) # this is >= 1
    A = Matrix.random(PolyRing, m, m, degree=deg)
    B = Matrix.random(PolyRing, m, n, degree=d-1)
    X = system_solve_expansion(A, B, d)
    BB = A*X  # should be B mod x^d
    matrix_truncate(BB, d)
    return B == BB

print "Testing expansion of linear system solution..."
correct = all([test_system_solve_expansion() for i in range(50)])
print " -->", "correct" if correct else "wrong"

def test_system_solve():
    p = next_prime(1000000, proof=True)  # using big prime to ensure Prob(A(0) singular) ~ 0
    PolyRing.<x> = GF(p)[]
    m = ZZ.random_element(3,10)
    deg = ZZ.random_element(20,40)
    d = m*deg
    A = Matrix.random(PolyRing, m, m, degree=deg)
    B = Matrix.random(PolyRing, m, 1, degree=d)
    X,f = system_solve(A, B)
    return A*X == f*B

print "Testing linear system solving..."
correct = all([test_system_solve() for i in range(50)])
print " -->", "correct" if correct else "wrong"


print "###########################"
print "# Testing MATRIX DIVISION #"
print "###########################"

def test_matrix_quo_rem():
    p = next_prime(1000000, proof=True)  # using big prime to ensure Prob(B row reduced) ~ 1
    PolyRing.<x> = GF(p)[]
    m = ZZ.random_element(3,10)
    n = ZZ.random_element(1,5)
    # A = random matrix of degree < d
    d = floor(m/n * 30) # this is >= 1
    A = Matrix.random(PolyRing, m, n, degree=d-1)
    # B = random matrix of row degree rdeg
    rdeg = [ZZ.random_element(20,40) for i in range(m)]
    B = Matrix(PolyRing, m, m)
    for i in range(m):
        B[i,:] = Matrix.random(PolyRing, 1, m, degree=rdeg[i])
    Q,R = matrix_quo_rem(A, B)
    if A != B*Q + R:
        return False
    rdegR = R.row_degrees()
    if any([rdegR[i] >= rdeg[i] for i in range(m)]):
        return False
    return True

print "Testing polynomial matrix division with remainder..."
correct = all([test_matrix_quo_rem() for i in range(50)])
print " -->", "correct" if correct else "wrong"


print "#############################"
print "# Testing APPROXIMANT BASIS #"
print "#############################"

def test_approximant_basis():
    PolyRing.<x> = GF(97)[]
    m = ZZ.random_element(1,6)
    n = ZZ.random_element(2,8)
    order = ZZ.random_element(5,100)
    shift = [ZZ.random_element(0,20) for k in range(n)]
    F = Matrix.random(PolyRing, m, n, degree=order-1)
    P,cdeg = approximant_basis(F, order, shift)
    # verify that cdeg is the shift-column degree
    if cdeg != P.column_degrees(shifts=shift):
        return False
    # verify that P is in shift-ordered weak Popov form
    if not P.is_weak_popov(row_wise=False, shifts=shift, ordered=True):
        return False
    # verify that F*P = 0 mod x^order
    residual = F*P
    certificate = copy(residual) # for not recomputing later
    matrix_truncate(residual, order)
    if residual != 0:
        return False
    # verify that the columns of P generate all approximants
    # --> cheap verification based on the form of the determinant of P and on a
    # constant matrix having full rank [see Giorgi-Neiger, ISSAC 2018]
    deg_det = sum([P[i,i].degree() for i in range(n)])
    if P.determinant() != (P(1).determinant() * x**deg_det):
        return False
    matrix_shift(certificate, -order)
    mat = Matrix.block([[certificate(0)], [P(0)]])
    if mat.rank() < n:
        return False
    return True

print "Testing minimal approximant basis..."
correct = all([test_approximant_basis() for i in range(15)])
print " -->", "correct" if correct else "wrong"

