#!/usr/bin/env sage


# As input we are given

# Target lattice
#L = {(u_{1}, ..., u_{n}) \in \zz^{n} |  \forall 1 \leq j \leq k: \prod_{i=1}^{n}(x - \alpha_{i})^{u_{i}} \equiv 1 \pmod{c_{j}(x)}}


# The following function takes as input c_{j}(x) and a list of \alpha_{i} and
# produces the representation of L stated above
# More specifically it gives a list of log_{\beta_j}(x- \alpha_{i}) j = 1, ..., n

# The algorithm proceeds as follows:
# 1.Find a generator of multiplicative groups (F_q[x] / c_j(x))*
# 2.Compute logs of x - \alpha_i with respect to it using any algorithm
# 3.output the list of log_{\beta_j}(x - \alpha_i) j = 1, ..., n


def parity_check_representation(field_params, modulus, alphas):
    #print "in PARITY CHECK"
    f = GF(field_params[0]**field_params[1], 'x', modulus = modulus)
    #The way we constructed it x is a generator

    logs = []
    for p in primes:
        log = compute_log(order, gen, p)
        logs.append(log)
    return vector(logs)


def compute_log(order, generator, element):
    #print "in LOG"
    return discrete_log(element, generator, order)


def dual_generating_set(modulus_factors, primes):
    #print "in DUAL GEN SET"
    dual_gens = identity_matrix(QQ, len(primes))
    for q, e in modulus_factors:
        euler_phi = q**e - q**(e-1)
        dual_gens = dual_gens.stack(parity_check_representation((q, e), primes)
                                    / euler_phi)
    return dual_gens


# input type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'
# input type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'
def hermite_form(B, indep_length):
    #print "in HERMITE FORM"
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
    basis = basis.matrix_from_rows(range(rows_number - indep_length, rows_number))
    return basis


# input type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'
# input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def primal_basis_from_dual(B):
    #print "in PRIMAL FROM DUAL"
    try:
        return (B * ~(B.T * B)).change_ring(ZZ)
    except Exception as msg:
        print "exception: ", msg


def test_lattice_construction(n, t):
    primes = primes_first_n(n+t)
    modulus_factors = primes[1:t+1]
    primes = [2] + primes[t+1:]
    modulus_factors = [(q, n) for q in modulus_factors]
    return (modulus_factors, primes, lattice_construction(modulus_factors, primes))


def lattice_construction(modulus_factors, primes):
    #time
    dual_gens = dual_generating_set(modulus_factors, primes)
    #time
    #dual_basis = hermite_form(dual_gens, len(primes))
    dual_basis = lll_wrap(dual_gens, len(primes))
    print dual_basis
    # matrix dimention is (n+t) * n
    #time
    primal_basis = primal_basis_from_dual(dual_basis)
    #determinant check
    phi_of_modulus = 1
    for (q, n) in modulus_factors:
        phi_of_modulus *= q**n - q**(n-1)
    #print abs(det(primal_basis)) == phi_of_modulus
    # In all test runs the result was an integer matrix,
    # no function threw an exception

    # rows are basis vectors here!! That's why .T
    return primal_basis.T


#print test_lattice_construction(18, 4)
