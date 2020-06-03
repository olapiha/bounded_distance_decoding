#!/usr/bin/env sage
import numpy as np
from sage.modules.free_module_integer import IntegerLattice

# As input we are given a lattice basis,
# prime numbers and modulus m used to construct it
# and a point point_t = (t_1, ..., t_n) in R^n such that point_t = x + e,
# where x - lattice point and e - bounded error

# The goal is to recover x from t.
# To achieve it we use algebraic sturucture of the lattice


#ToDo NEED TO FIND A WAY TO WORK WITH POLYNOMIALS MODULO C 

def positive_discrete_error(alphas, modulus, Fx, point_t):
    print "in POSITIVE_DISCRETE_ERROR"
    n = len(alphas)
    #Fxmod = Fx.quotient(modulus)
    prod_modulo = Fx(1)
    #gen = Fxmod(x)
    gen = Fx(x)
    print "making primes"
    primes = [gen - alpha for alpha in alphas]
    for i in range(n):
        print "in first for loop"
        prime_power = Fx(1)
        for j in range(point_t[i]):
            print prime_power * modulus.quo_rem(primes[i])[1]
            print modulus.quo_rem(prime_power * modulus.quo_rem(primes[i])[1])
            prime_power = modulus.quo_rem(prime_power * modulus.quo_rem(primes[i])[1])[1]
        prod_modulo = modulus.quo_rem(prod_modulo * prime_power)[1]
    print "starting factorization"
    error = []
    for p in primes:
        e = 0
        while prod_modulo.quo_rem(p)[1] == 0:
            prod_modulo = prod_modulo.quo_rem(p)[0]
            e += 1
        error.append(e)
    error = vector(error)
    return error


def discrete_error(alphas, modulus, Fx, point_t):
    n = len(alphas)
    Fxmod = Fx.quotient(modulus)
    prod_modulo = Fxmod(1)
    gen = Fxmod(x)
    primes = [gen - alpha for alpha in alphas]
    for i in range(n):
        prod_modulo = Fxmod(prod_modulo * Fxmod(primes[i])**point_t[i])
    (numerator, denominator) = rational_function_reconstruction(prod_modulo, modulus, Fx)
    #print numerator.factor()
    #print denominator.factor()
    num_error = []
    for p in primes:
        e = 0
        while prod_modulo.quo_rem(p)[1] == 0:
            prod_modulo = prod_modulo.quo_rem(p)[0]
            e += 1
        num_error.append(e)
    num_error = vector(num_error)
    den_error = []
    for p in primes:
        e = 0
        while prod_modulo.quo_rem(p)[1] == 0:
            prod_modulo = prod_modulo.quo_rem(p)[0]
            e += 1
        den_error.append(e)
    den_error = vector(den_error)
    overall_error = num_error - den_error
    return overall_error


def rational_function_reconstruction(g, f, Fx):
    assert f.degree() > g.degree()
    if g == 0:
        assert f.degree() > 3 # do I actually need this? Maybe just for coherence
        return (0,1)
    (r0, r1) = (f, g)
    (t0, t1) = (0, 1)
    q = Fx(1)
    while q.degree() <= f.degree()/2:
        (q, r) = r0.quo_rem(r1)
        (r0, r1) = (r1, r) # r = r0 - q*r1
        (t0, t1) = (t1, t0 - q*t1)
    assert (r0.degree() + t0.degree() < f.degree()/2)
    assert r0.gcd(t0)
    return (r0, t0)


def find_generator(ring):
    #print "in FIND_GENERATOR"
    while True:
        candidate = ring.random_element()
        desired_order = ring.order()-1
        factors = list(factor(desired_order))
        check = True
        for (fac, foo) in factors:
            possible_order = desired_order // fac
            if candidate**possible_order == 1:
                check = False
                break
        if check:
            break
    return candidate


# B must be an integer
def test_decoding(q, d, n, k, B):
    load("lattice computation polynomials.sage")
    (alphas, modulai, Fx, basis) = test_lattice_construction(q, d, n, k)
    modulus = Fx(1)
    for c in modulai:
        modulus *= c

    # generate coordinates of a lattice point in the basis above
    coordinates = vector(np.random.randint(0, 10, n))
    lattice_point = basis * coordinates

    # generate integer noise of l_1 norm <= B
    while True:
        noise = vector(np.random.randint(0, B+1, n))
        if (noise.norm(1) <= B): break
    print "norm of the noise: ", noise.norm()
    point_t = lattice_point + noise
    return (positive_discrete_error(alphas, modulus, Fx, point_t), noise)



result = test_decoding(3 ^ 5, 2, 10, 6, 2)
print result[0]
print result[1]
