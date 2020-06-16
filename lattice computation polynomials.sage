#!/usr/bin/env sage
import numpy as np

# As input we are given

# Target lattice
# L = {(u_{1}, ..., u_{n}) \in \zz^{n} |  \forall 1 \leq j \leq k: \prod_{i=1}^{n}(x - \alpha_{i})^{u_{i}} \equiv 1 \pmod{c_{j}(x)}}


# The following function takes as input c_{j}(x) and a list of \alpha_{i} and
# produces the representation of L stated above
# More specifically it gives a list of log_{\beta_j}(x- \alpha_{i}) i = 1, ..., n

# The algorithm proceeds as follows:
# 1.Find a generator of multiplicative groups (F_q[x] / c_j(x))*
# 2.Compute logs of x - \alpha_i with respect to it using any algorithm
# 3.output the list of log_{\beta_j}(x - \alpha_i) j = 1, ..., n

#Note: Basis vectors are rows until the output in the very end.

def discrete_logs_vector(field_order, modulus, alphas):
    #print("in DISCRETE LOGS VECTOR")
    F.<y> = GF(field_order)
    Fx.<x> = PolynomialRing(F)
    I = modulus * Fx
    Fx_quotient = Fx.quotient(I)
    generator = find_generator(Fx_quotient)
    #print("out FIND_GENERATOR")
    logs = []
    for alpha in alphas:
        element = x - alpha
        log = discrete_log(element,  generator, Fx_quotient.order()-1)
        logs.append(log)
    return vector(logs)


def find_generator(ring):
    #print("in FIND_GENERATOR")
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

# Here we produce a generating set of row vectors.
def dual_generating_set(field_order, modulai, alphas):
    #print("in DUAL GEN SET")
    d = modulai[0].degree()
    dual_gens = identity_matrix(QQ, len(alphas))
    for modulus in modulai:
        dual_gens = dual_gens.stack(discrete_logs_vector(
            field_order, modulus, alphas) / (field_order^d - 1))
    return dual_gens


def hermite_form(B, indep_length):
    #print("in HERMITE FORM")
    den = B.denominator()
    int_B = (den*B).change_ring(ZZ)
    (H, U) = int_B.hermite_form(transformation=True)
    basis = (1/den)*H
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows last and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(range(indep_length))
    return basis

# LLL takes as input a set of row vectors
def lll_wrap(B, indep_length):
    #print("in LLL")
    basis = B.LLL()
    rows_number = len(B.rows())
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows first and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(
        range(rows_number - indep_length, rows_number))
    return basis


def primal_basis_from_dual(B):
    #print("in PRIMAL FROM DUAL")
    return (B * ~(B.T * B)).change_ring(ZZ)


def parity_check_representation(field_order, modulai, alphas):
    #print("in PARITY CHECK")
    parity_check = matrix(discrete_logs_vector(field_order, modulai[0], alphas))
    for modulus in modulai[1:]:
        parity_check = parity_check.stack(discrete_logs_vector(field_order, modulus, alphas))
    return parity_check


# Don't forget: q should be prime or prime power!
# n < q
# k < q^d / d
# q field order
# d degree of polynomials
# n dimention of the lattice
# k number of irreducible c_j(x)
def check_parameters(q, d, n, k):
    assert len(q.factor()) == 1
    assert n <= q
    assert k < q^d /d


def test_lattice_construction(field_order, d, n, k):
    check_parameters(field_order, d, n, k)
    F.<y> = GF(field_order)
    Fx.<x> = PolynomialRing(F)
    modulai = Set([])
    while modulai.cardinality() < k:
        modulai += Set([Fx.irreducible_element(d, algorithm="random")])
    alphas = Set([])
    while alphas.cardinality() < n:
        alphas = alphas + Set([F.random_element()])
    basis = lattice_construction_v2(field_order, list(modulai), list(alphas))
    basis1 = lattice_construction(field_order, list(modulai), list(alphas), None)
    return (alphas, modulai, basis, basis1)


def lattice_construction(field_order, modulai, alphas, parity_check):
    if parity_check is None:
        dual_gens = dual_generating_set(field_order, modulai, alphas)
    else:
        d = modulai[0].degree()
        dual_gens = identity_matrix(QQ, len(alphas)).stack(parity_check / (field_order^d - 1))
    dual_basis = hermite_form(dual_gens, len(alphas))
    #dual_basis = lll_wrap(dual_gens, len(alphas))
    # matrix dimention is (n+t) * n
    primal_basis = primal_basis_from_dual(dual_basis)
    # determinant check
    det_upperbound = (field_order ^ (modulai[0].degree()) - 1) ^ len(modulai)
    #print("determinant check:")
    #print("abs value of det(basis): ", abs(primal_basis.det()))
    #print("upperbound:              ", det_upperbound)
    # rows are basis vectors here but want columns!! That's why .T
    return primal_basis.T

def lattice_construction_v2(field_order, modulai, alphas):
    parity_check = matrix(QQ, parity_check_representation(field_order, modulai, alphas))
    echelon_form_in_QQ =  parity_check.echelon_form()
    d = modulai[0].degree()
    try:
        order = field_order^d -1
        systematic_form = matrix(Zmod(order), echelon_form_in_QQ)
        k = len(modulai)
        new_columns = [(-i, -j) for (i, j) in systematic_form.columns()[k:]] + systematic_form.columns()[:k]
        generator_mod_q = matrix(QQ, new_columns).T
        generator = (identity_matrix(QQ, len(alphas))*order).stack(generator_mod_q)
        return hermite_form(generator, len(alphas))
    except Exception as e:
        print("Falied to reduce to systematic form")
        print("Running the slower algorithm")
        return lattice_construction(field_order, list(modulai), list(alphas), parity_check)




res = test_lattice_construction(7, 2, 3, 2)
#print res[0]
print res[2]
print res[3]
