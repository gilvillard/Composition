reset()
load('utils.sage')
load('balanced_basis.sage')
load('inverse_composition.sage')
load('modular_composition.sage')

PolyRing.<y> = GF(997)[]
PolyRingX.<x> = GF(997)[]

n=2000
m=ceil(n^0.2)
# we choose exponent 0.2 instead of 1/omega (~0.33 in Sage), since
# the used prototype implementation of reconstruction of recurrent matrix
# sequence is not very efficient

# Balanced basis
g = PolyRing.random_element(degree=n).monic()
a = PolyRing.random_element(degree=n-1)
(B, pow_a) = balanced_basis(g, a, m, store_powers=True, verbose=True)
# with store_powers=True, it also returns pow_a which contains a^k mod g for k in range(2ceil(n/m))

# define the quotient ring K[y] / <g>
QuoRing.<yy> = PolyRing.quotient(g)

# Inverse composition
b = PolyRing.random_element(degree=n-1)
h = inverse_composition(g, a, b, m, verbose=True)

if h(a(yy)) == b(yy):
    print "Inverse composition --> Correct"
else:
    print "Inverse composition --> Wrong"


# Modular composition
#h = PolyRingX.random_element(degree=n-1)
#b = modular_composition(g, h, a, m, verbose=True)
