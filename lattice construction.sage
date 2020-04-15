#!/usr/bin/env sage

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


def parity_check_representation((q, e), primes):
    #print "in PARITY CHECK"
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
    #print "in LOG"
    return discrete_log(element, generator, order)


def dual_generating_set(modulus_factors, primes):
    print "in DUAL GEN SET"
    dual_gens = identity_matrix(QQ, len(primes))
    for q, e in modulus_factors:
        euler_phi = q**e - q**(e-1)
        dual_gens = dual_gens.stack(parity_check_representation((q, e), primes)
                                    / euler_phi)
    return dual_gens


# input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
# output type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def lll_wrapper(B):
    print "in LLL"
    return B.LLL()


# input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
# input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def primal_basis_from_dual(B):
    print "in PRIMAL FROM DUAL"
    try:
        return (B * ~(B.T * B)).change_ring(ZZ)
    except Exception as msg:
        print "exception: ", msg


def test_main(n, t):
    primes = primes_first_n(n+t)
    modulus_factors = primes[1:t+1]
    primes = [2] + primes[t+1:]
    modulus_factors = [(q, n) for q in modulus_factors]
    print main(modulus_factors, primes)


def main(modulus_factors, primes):
    time dual_gens = dual_generating_set(modulus_factors, primes)
    time dual_basis = lll_wrapper(dual_gens)
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # lll outputs zero rows first and we know there will be exactly t of them so:
    dual_basis = dual_basis.delete_rows(range(len(modulus_factors)))
    # matrix dimention is (n+t) * n
    time primal_basis = primal_basis_from_dual(dual_basis)
    # In all test runs the result was an integer matrix,
    # no function threw an exception
    return primal_basis


test_main(10, 4)
