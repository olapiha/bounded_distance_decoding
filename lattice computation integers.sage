#!/usr/bin/env sage
import time


# As input we are given an interger m and its factorization
# m = \Pi_{1}^{t} q_{i}^{e_{i}}, q_i its prime factors
# We are also given a set of p_{i} 1 \leq i \leq n ->
# prime numbers that don't divide m and bounded by some constant B

# Target lattice L a is an intersection of
# L_i = ((x_1, ..., x_n) \in Z | \sum_{1}^{n} x_j * log_{\beta_i}(p_j) = 0 \pmod(\phi(q_{i}^{e_i})))
# The following function takes as input q_{i}^{e_i} and
# produces the representation of L_i stated above
# More specifically it gives a list of log_{\beta_i}(p_j) j = 1, ..., n

# The algorithm proceeds as follows:
# 1.Find a generator of multiplicative groups Z / q_{i}^{e_i}}Z
# 2.Compute logs of p_j with respect to it using PH & Pollard-rho algorithms
# 3.output the list of log_{\beta_i}(p_j) j = 1, ..., n


def parity_check_representation(q, e, primes):
    order = q**e - q**(e-1)
    q = q**e
    Zq = Zmod(q)
    multiplicative_group = Zq.unit_group()
    # I checked on few examples the group is actually cyclic
    gen = multiplicative_group.gens_values()[0]
    logs = []
    for p in primes:
        log = compute_log(order, gen, p)
        logs.append(log)
    return vector(logs)


def compute_log(order, generator, element):
    return discrete_log(element, generator, order)


def dual_generating_set(modulus_factors, primes):
    dual_gens = identity_matrix(QQ, len(primes))
    for q, e in modulus_factors:
        euler_phi = q**e - q**(e-1)
        dual_gens = dual_gens.stack(parity_check_representation(q, e, primes)
                                    / euler_phi)
    return dual_gens


def hermite_form(B, indep_length):
    den = B.denominator()
    int_B = (den*B).change_ring(ZZ)
    (H, U) = int_B.hermite_form(transformation = True)
    basis = (1/den)*H
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows last and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(range(indep_length))
    return basis


def lll_wrap(B, indep_length):
    basis = B.LLL()
    rows_number = len(B.rows())
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows first and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(
        range(rows_number - indep_length, rows_number))
    return basis


def primal_basis_from_dual(B):
    return (B * ~(B.T * B)).change_ring(ZZ)


def test_lattice_construction(n, t):
    primes = primes_first_n(n+t)
    modulus_factors = primes[n:]
    primes = primes[:n]
    modulus_factors = [(q, round((ln(3)* n)/(ln(q) * t))) for q in modulus_factors]
    return (modulus_factors, primes, lattice_construction(modulus_factors, primes))


def lattice_construction(modulus_factors, primes):
    begin = time.time()
    dual_gens = dual_generating_set(modulus_factors, primes)
    print(time.time() - begin)
    #dual_basis = hermite_form(dual_gens, len(primes))
    begin = time.time()
    dual_basis = lll_wrap(dual_gens, len(primes))
    print(time.time() - begin)
    # matrix dimention is (n+t) * n
    begin = time.time()
    primal_basis = primal_basis_from_dual(dual_basis)
    print(time.time() - begin)
    #determinant check
    phi_of_modulus = 1
    for (q, n) in modulus_factors:
        phi_of_modulus *= q**n - q**(n-1)
    #print(abs(det(primal_basis)) == phi_of_modulus)
    # In all test runs the result was an integer matrix,
    # no function threw an exception

    # rows are basis vectors here!! That's why .T
    return primal_basis.T

'''
begin = time.time()
print(test_lattice_construction(200, 20))
print(time.time() - begin)
'''
