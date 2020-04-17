#!/usr/bin/env sage
import numpy
from sage.modules.free_module_integer import IntegerLattice

# As input we are given a lattice basis,
# prime numbers and modulus m used to construct it
# and a point point_t = (t_1, ..., t_n) in R^n such that point_t = x + e,
# where x - lattice point and e - bounded error

# The goal is to recover x from t and test it with real lattice instance
# To achieve it we use algebraic sturucture of the lattice


def positive_discrete_error(primes, modulus, point_t):
    prod_modulo = 1
    n = len(primes)
    Zm = Zmod(modulus)
    for i in range(n):
        prod_modulo = Zm(prod_modulo * Zm(primes[i]**point_t[i]))
    error = []
    for p in primes:
        e = 0
        Zp = Zmod(p)
        while Zp(prod_modulo) == 0:
            prod_modulo = prod_modulo // p
            e += 1
        error.append(e)
    error = vector(error)
    return error


def discrete_error(primes, modulus, point_t):
        prod_modulo = 1
        n = len(primes)
        Zm = Zmod(modulus)
        for i in range(n):
            prod_modulo = Zm(prod_modulo * Zm(primes[i]**point_t[i]))

        (numerator, denominator) = rational_number_reconstruction(prod_modulo, modulus)
        print primes
        print numerator.factor()
        print denominator.factor()
        numerator_error = []
        for p in primes:
            e = 0
            Zp = Zmod(p)
            while Zp(numerator) == 0:
                numerator = numerator // p
                e += 1
            numerator_error.append(e)
        numeator_error = vector(numerator_error)

        denominator_error = []
        for p in primes:
            e = 0
            Zp = Zmod(p)
            while Zp(denominator) == 0:
                denominator = denominator // p
                e += 1
            denominator_error.append(e)
        denominator_error = vector(denominator_error)
        overall_error = numeator_error - denominator_error
        return overall_error


def rational_number_reconstruction(f, modulus):
    Zm = Zmod(modulus)
    basis = IntegerLattice([[Zm(f),1],[modulus,0]])
    return basis.shortest_vector()


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
    return (positive_discrete_error(primes, modulus, point_t), noise)


result = test_decoding(50, 3, 2)
print result[0]
print result[1]
