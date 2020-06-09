#!/usr/bin/env sage
import numpy as np

# As input we are given a lattice basis,
# prime numbers and modulus m used to construct it
# and a point point_t = (t_1, ..., t_n) in R^n such that point_t = x + e,
# where x - lattice point and e - bounded error

# The goal is to recover x from t.
# To achieve it we use algebraic sturucture of the lattice


def positive_discrete_error(alphas, modulus, Fx, point_t):
    n = len(alphas)
    Fxmod = Fx.quotient(modulus)
    prod_modulo = Fxmod(1)
    gen = Fxmod(x)
    primes = [gen - alpha for alpha in alphas]
    for i in range(n):
        prod_modulo = Fxmod(prod_modulo * Fxmod(primes[i])**point_t[i])
    error = []
    prod_modulo = lift(prod_modulo)
    for alpha in alphas:
        e = 0
        while prod_modulo(alpha) == 0:
            prod_modulo = prod_modulo.quo_rem(lift(gen-alpha))[0]
            e += 1
        error.append(e)
    error = vector(error)
    return error


def discrete_error(alphas, modulus, Fx, point_t):
    n = len(alphas)
    Fxmod = Fx.quotient(modulus)
    prod_modulo = Fxmod(1)
    gen = Fx(x)
    primes = [gen - alpha for alpha in alphas]
    for i in range(n):
        prod_modulo = Fxmod(prod_modulo * Fxmod(primes[i])**point_t[i])
    prod_modulo = lift(prod_modulo)
    (numerator, denominator) = rational_function_reconstruction(prod_modulo, modulus, Fx)
    num_error = []
    for alpha in alphas:
        e = 0
        while numerator(alpha) == 0:
            numerator = numerator.quo_rem(gen-alpha)[0]
            e += 1
        num_error.append(e)
    num_error = vector(num_error)
    den_error = []
    for alpha in alphas:
        e = 0
        while denominator(alpha) == 0:
            denominator = denominator.quo_rem(gen-alpha)[0]
            e += 1
        den_error.append(e)
    den_error = vector(den_error)
    overall_error = num_error - den_error
    return overall_error


def rational_function_reconstruction(g, f, Fx):
    assert f.degree() > g.degree()
    if g == 0:
        assert f.degree() > 3 # for coherence
        return (0,1)
    (r0, r1) = (f, g)
    (t0, t1) = (Fx(0), Fx(1))
    q = Fx(1)
    while q.degree() <= f.degree()/2:
        (q, r) = r0.quo_rem(r1)
        (r0, r1) = (r1, r) # r = r0 - q*r1
        (t0, t1) = (t1, t0 - q*t1)
    assert (r0.degree() + t0.degree() < f.degree()/2)
    assert r0.gcd(t0)
    return (r0, t0)


# B must be an integer
def test_decoding(q, d, n, k, B):
    load("lattice computation polynomials.sage")
    (alphas, modulai, basis) = test_lattice_construction(q, d, n, k)
    F.<y> = GF(q)
    Fx.<x> = PolynomialRing(F)
    modulus = Fx(1)
    for c in modulai:
        modulus *= c
    # generate coordinates of a lattice point in the basis above
    coordinates = vector(np.random.randint(0, 10, n))
    lattice_point = basis * coordinates
    # generate integer noise of l_1 norm <= B
    while True:
        noise = vector(np.random.randint(-B, B+1, n))
        if (noise.norm(1) <= B): break
    print "l1 norm of the noise: ", noise.norm(1)
    point_t = lattice_point + noise
    return (discrete_error(alphas, modulus, Fx, point_t), noise)


for i in range(100):
    result = test_decoding(3 ^ 7, 2, 7, 6, 2)
    print "calculated error: ", result[0]
    print "real error:       ", result[1]
