#!/usr/bin/env sage

#As input we are given an interger m and its factorization m = \Pi_{1}^{t} q_{i}^{e_{i}}, q_i its prime factors
#We are also given a set of p_{i} 1 \leq i \leq n ->
#prime numbers that don't divide m and bounded by some constant B

#Target lattice L a is an intersection of
#L_i = ((x_1, ..., x_n) \in Z | \sum_{1}^{n} x_j * log_{\beta_i}(p_j) = 0 \pmod(\phi(q_{i}^{e_i})))
#The following function takes as input q_{i}^{e_i} and produces the representation of L_i stated above
#More specifically it gives a list of log_{\beta_i}(p_j) j = 1, ..., n

#The algorithm proceeds as follows:
#1.Find a generator of multiplicative groups Z / q_{i}^{e_i}}Z
#2.Compute logs of p_j with respect to it using PH & Pollard-rho algorithms
#3.output the list of log_{\beta_i}(p_j) j = 1, ..., n
def parity_check_representation((q, e), primes):
    order = q**e - q**(e-1)
    q = q**e
    Zq = Zmod(q)
    multiplicative_group = Zq.unit_group()
    #Q: Should I check if it is actually cyclic?
    gen = multiplicative_group.gens_values()[0]
    print gen
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
        dual_gens = dual_gens.stack(parity_check_representation((q,e), primes) / (q**e - q**(e-1)))
    return dual_gens


#input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
#output type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def lll_wrapper(B):
    return B.LLL()


#input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
#input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def primal_basis_from_dual(B):
    try:
        return (B * ~(B.T * B)).change_ring(ZZ)
    except Exception as msg:
        print "exception: ", msg


def test_pbfd(n):
    rows = []
    for i in range(n):
        print "Input", i, "th row here"
        ith_row = []
        for j in range(n):
            element = int(input("next element:"))
            ith_row.append(element)
        rows.append(ith_row)
    B = matrix(ZZ, rows)
    print "dual basis:"
    print B
    print "primal basis: "
    print primal_basis_from_dual(B)


def test_lll(n):
    rows = []
    for i in range(n):
        print "Input", i, "th row here"
        ith_row = []
        for j in range(n):
            element = int(input("next element:"))
            ith_row.append(element)
        rows.append(ith_row)
    B = matrix(ZZ, rows)
    print "initial basis:"
    print B
    print "reduced basis: "
    print lll_wrapper(B)


def test_pcr(modulus_factors, primes):
    res = []
    for factor in modulus_factors:
        res.append(parity_check_representation(factor, primes))
    print res


def test_main(n, t):
    primes = primes_first_n(n)
    lowerbound = primes[n-1]
    modulus_factors = []
    for i in range(t):
        q = random_prime(1000, proof=None, lbound=lowerbound+1)
        e = ZZ.random_element(1,n)
        modulus_factors.append((q, e))

    print modulus_factors
    print primes
    print main(modulus_factors, primes)


def main(modulus_factors, primes):
    dual_gens = dual_generating_set(modulus_factors, primes)
    dual_basis = lll_wrapper(dual_gens)
    primal_basis = primal_basis_from_dual(dual_basis)
    primal_basis = primal_basis.delete_rows(range(len(modulus_factors)))
    return primal_basis


test_main(50, 4)
