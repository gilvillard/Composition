reset()
load('utils.sage')
load('balanced_basis.sage')
load('inverse_composition.sage')
load('modular_composition.sage')
load('power_projections.sage')

PolyRing.<y> = GF(997)[]
PolyRingX.<x> = GF(997)[]

n=400
m=ceil(n^0.2)
# we choose exponent 0.2 instead of 1/omega (~0.33 in Sage), since
# the used prototype implementation of reconstruction of recurrent matrix
# sequence is not very efficient

# BALANCED BASIS
g = PolyRing.random_element(degree=n).monic()
a = PolyRing.random_element(degree=n-1)
(B, pow_a) = balanced_basis(g, a, m, store_powers=True, verbose=True)
# with store_powers=True, it also returns pow_a which contains a^k mod g for k in range(2ceil(n/m))

# define the quotient ring K[y] / <g>
QuoRing.<yy> = PolyRing.quotient(g)

# INVERSE COMPOSITION
b = PolyRing.random_element(degree=n-1)
h = inverse_composition(g, a, b, m, verbose=True)

if h(a(yy)) == b(yy):
    print "Inverse composition --> Correct"
else:
    print "Inverse composition --> Wrong"

# MODULAR COMPOSITION
h = PolyRingX.random_element(degree=n-1)
b = modular_composition(g, a, h, m, verbose=True)

if h(a(yy)) == b(yy):
    print "Composition --> Correct"
else:
    print "Composition --> Wrong"

# POWER PROJECTIONS
elly = [GF(997).random_element() for k in range(n)]
ella = power_projections(g, a, elly, m, verbose=True)

print "naive recomputation of power projections for testing purposes..."
correct = (ella[0] == elly[0])
k=1
ak = a(yy)
while correct and k<n:
    # list of coefficients, padded with zeroes if necessary
    coeffs = list(ak)
    coeffs = coeffs + [0]*(n-len(coeffs))
    ellak = (Matrix(GF(997),1,n,elly) * Matrix(GF(997),n,1,coeffs))[0,0]
    correct = (ella[k] == ellak)
    ak = ak * a(yy)
    k = k+1

if correct:
    print "PowerProjections --> Correct"
else:
    print "PowerProjections --> Wrong"
