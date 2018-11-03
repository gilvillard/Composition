reset()
load('utils.sage')
load('balanced_basis.sage')
load('inverse_composition.sage')
load('modular_composition.sage')

PolyRing.<y> = GF(997)[]
PolyRingX.<x> = GF(997)[]

n=110
m=ceil(n^0.33)

# Balanced basis
#g = PolyRing.random_element(degree=n).monic()
#a = PolyRing.random_element(degree=n-1)
#(B, pow_a) = balanced_basis(g, a, m, store_powers=True, verbose=True)

# Inverse composition
#g = PolyRing.random_element(degree=n).monic()
#a = PolyRing.random_element(degree=n-1)
#b = PolyRing.random_element(degree=n-1)
n = 4
m = 2
g = y^4
a = y^2
b = PolyRing(y)
(h, B) = inverse_composition(g, a, b, m, verbose=True)
if h(a) % g == b:
    print "Inverse composition: correct"
else:
    print "Inverse composition: error"


#if False:
#    # Check same ideal (incidentally, tests gg == B.det() == charpoly(a mod g))
#    H = hermite_form(B)
#    gg = H[0,0]
#    aa = -H[0,1]
#    bivPolyRing.<X,Y> = GF(997)[]
#    I1 = bivPolyRing.ideal(g(Y), X-a(Y))
#    I2 = bivPolyRing.ideal(gg(X), Y-aa(X))
#    print "Change of order lex(x<y) --> lex(x>y) ?", I1 == I2
#
#    # Inverse composition
#    b = PolyRing.random_element(degree=n-1)
#    h = b(aa(x)) % gg(x)
#    print "Inverse composition via composition ?", b == h(a) % g

## Modular composition
#h = PolyRingX.random_element(degree=n-1)
#b = modular_composition(g, h, a, m, verbose=True)
