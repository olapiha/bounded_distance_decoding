#!/usr/bin/env sage
import numpy as np
import time
import random


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
    print("in DISCRET_ERROR")
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
    print("in RFR")
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
    start = time.time()
    (alphas, modulai, basis) = test_lattice_construction(q, d, n, k)
    end = time.time()
    print("lattice construction time:", end-start)
    F.<y> = GF(q)
    Fx.<x> = PolynomialRing(F)
    modulus = Fx(1)
    for c in modulai:
        modulus *= c
    # generate coordinates of a lattice point in the basis above
    coordinates = vector(np.random.randint(0, 10, n))
    lattice_point = basis * coordinates
    # generate integer noise of l_1 norm <= B
    non_zero_indices = list(np.random.choice(range(n), B))
    while True:
        noise = []
        for i in range(n):
            if i in non_zero_indices:
                noise.append(random.randint(-1, 2))
            else:
                noise.append(0)
        noise = vector(noise)
        if (noise.norm(1) <= B) and (0 < noise.norm(1)): break
    '''
    noise = vector([1 for i in range(B//2)] + [-1 for i in range(B//2)] + [0 for i in range(n-48)])
    '''
    print("l1 norm of the noise: ", noise.norm(1))
    point_t = lattice_point + noise
    start = time.time()
    e = discrete_error(alphas, modulus, Fx, point_t)
    end = time.time()
    print("decoding time:", end-start)
    return (e, noise)


def test_decoding_with_optimal_parameters(n):
    q = n.next_prime()
    d = 2
    k = (n/(2*ln(n))).round()
    B = (d*k/(2*(q^d-1)^(k/n))).round()
    print(q,d,k,B)
    return test_decoding(q, d, n, k, B)

'''
start = time.time()
test_decoding_with_optimal_parameters(50)
end = time.time()
print("execution time": end-start)
'''
