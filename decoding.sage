#!/usr/bin/env sage
import numpy

# As input we are given a lattice basis,
# prime numbers and modulus m used to construct it
# and a point point_t = (t_1, ..., t_n) in R^n such that point_t = x + e,
# where x - lattice point and e - bounded error

# The goal is to recover x from t and test it with real lattice instance
# To achieve it we use algebraic sturucture of the lattice


def positive_discrete_error(primes, modulus, point_t):
    prod_modulo = 1
    n = len(primes)
    R = Zmod(modulus)
    for i in range(n):
        prod_modulo = R(prod_modulo * R(primes[i]**point_t[i]))
    error = []
    for p in primes:
        e = 0
        R = Zmod(p)
        while R(prod_modulo) == 0:
            prod_modulo = prod_modulo // p
            e += 1
        error.append(e)
    error = vector(error)
    return point_t - error


# B must be an integer
def test_decoding(n, t, B):
    load("lattice construction.sage")
    (modulus_factors, primes, basis) = test_lattice_construction(n, t)
    modulus = 1
    for (q, e) in modulus_factors:
        modulus *= q**e

    # generate coordinates of a lattice point in the basis above
    coordinates = vector(numpy.random.randint(0, 10, n))
    lattice_point = basis * coordinates

    # generate noise coords in default basis of Z^n
    noise = vector(numpy.random.randint(0, B, n))
    point_t = lattice_point + noise
    return (positive_discrete_error(primes, modulus, point_t), lattice_point)


#result = test_decoding(10, 3, 2)
#print result[0]
#print result[1]
