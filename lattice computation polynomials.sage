#!/usr/bin/env sage


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


def parity_check_representation(field_order, modulus, alphas):
    #print "in PARITY CHECK"
    #print generator
    F.<y> = GF(field_order)
    Fx.<x> = PolynomialRing(F)
    I = modulus * Fx
    Fx_quotient = Fx.quotient(I)
    xbar = Fx_quotient(x)
    generator = find_generator(Fx_quotient)
    logs = []
    for alpha in alphas:
        element = xbar - alpha
        log = discrete_log(element,  generator, Fx_quotient.order()-1)
        logs.append(log)
    return vector(logs)


def find_order(element):
    piv = element ^ 2
    order = 2
    while True:
        if piv == element: break
        piv = piv * element
        order += 1
    return order-1

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


def dual_generating_set(field_order, modulai, alphas):
    #print "in DUAL GEN SET"
    d = modulai[0].degree()
    dual_gens = identity_matrix(QQ, len(alphas))
    for modulus in modulai:
        dual_gens = dual_gens.stack(parity_check_representation(
            field_order, modulus, alphas) / (field_order^d - 1))
    return dual_gens


def hermite_form(B, indep_length):
    #print "in HERMITE FORM"
    den = B.denominator()
    int_B = (den*B).change_ring(ZZ)
    (H, U) = int_B.hermite_form(transformation=True)
    basis = (1/den)*H
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows last and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(range(indep_length))
    return basis


def lll_wrap(B, indep_length):
    #print "in LLL"
    basis = B.LLL()
    rows_number = len(B.rows())
    # Let  n= number of primes, t= number of factors
    # number of generators is n+t but the dimention of the lattice is n
    # hermite_form outputs zero rows first and we know there will be exactly t of them so:
    basis = basis.matrix_from_rows(
        range(rows_number - indep_length, rows_number))
    return basis


# input type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'
# input type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'
def primal_basis_from_dual(B):
    #print "in PRIMAL FROM DUAL"
    try:
        return (B * ~(B.T * B)).change_ring(ZZ)
    except Exception as msg:
        print "exception: ", msg

# Don't forget: q should be prime or prime power!
# n < q
# k < q^d / d
# q field order
# d degree of polynomials
# n dimention of the lattice
# k number of irreducible c_j(x)
def check_parameters(q, d, n, k):
    assert len(q.factor()) == 1
    assert n < q
    assert k < q^d /d


def test_lattice_construction(q, d, n, k):
    check_parameters(q, d, n, k)
    field_order = q
    F.<y> = GF(field_order)
    Fx.<x> = PolynomialRing(F)
    modulai = Set([])
    while True:
        modulai += Set([Fx.irreducible_element(d, algorithm="random")])
        if modulai.cardinality() == k:
            break
    alphas = Set([])
    while True:
        alphas = alphas + Set([F.random_element()])
        if alphas.cardinality() == n:
            break
    alphas_not_working = set([y^3 + y + 2, y^3 + 2*y, 2*y^4 + 2*y^3 + y^2 + 2*y + 2, 2*y^4 + y^2 + 1, y^2 + y + 2, 2*y^4 + y^3 + 2*y^2, 2*y^4 + 2*y^3 + 2*y^2 + 2*y, 2*y^4 + y^3 + 2*y^2 + y + 1, 2*y^2 + 2*y + 1, y^3 + 2*y^2 + 2])
    modulai_not_working = set([x^2 + (y^2 + 1)*x + 2*y^3 + y^2 + y, x^2 + (y^3 + 2*y^2 + 2*y + 2)*x + 2*y^3 + 1, x^2 + (2*y^4 + 2*y^3 + 1)*x + y^4 + y^3 + 2*y^2 + 2*y + 2, x^2 + (2*y^4 + y^2 + 2)*x + 2*y^3 + y^2 + y + 1, x^2 + (2*y^3 + y^2 + 2*y + 1)*x + 2*y^4 + y^3 + 2*y, x^2 + (y^4 + 2*y^3 + 2*y^2 + 2*y)*x + 2*y^2 + y + 1])
    alphnw = set([2*y^3 + y^2 + 2, 2*y^4, 1, y^3 + y^2 + 2*y + 2, 2*y^4 + 2*y^3 + 2*y^2, 2*y^4 + y^2 + 1, 2*y^4 + 2*y^3 + 2*y^2 + 2, 2*y^3 + 1, y^4 + 2*y^3 + y, y^4 + y^2 + y + 2])
    modnw = set([x^2 + (2*y^4 + y^3)*x + 2*y^3 + 2*y^2 + 1, x^2 + (2*y^3 + 2*y^2 + 2*y)*x + y^4 + 2*y^3 + 2*y^2 + y + 1, x^2 + (y^3 + 2*y^2 + 2*y + 2)*x + y^3 + y^2 + 1, x^2 + (y^4 + 2*y^3 + y^2 + 2*y + 2)*x + 2*y^4 + 2*y^3 + y^2 + 2*y, x^2 + (2*y^3 + y)*x + 2*y^3 + y^2 + 2*y + 2, x^2 + (y^4 + y + 2)*x + 2*y^3 + 2*y^2 + 2*y + 2])
    return  (alphas, modulai, lattice_construction(field_order, list(modnw), list(alphnw)))


def lattice_construction(field_order, modulai, alphas):
    dual_gens = dual_generating_set(field_order, modulai, alphas)
    #dual_basis = hermite_form(dual_gens, len(alphas))
    dual_basis = lll_wrap(dual_gens, len(alphas))
    # matrix dimention is (n+t) * n
    primal_basis = primal_basis_from_dual(dual_basis)
    # determinant check
    supposed_determinant = (field_order ^ (
        modulai[0].degree()) - 1) ^ len(modulai)
    print "determinant check:"
    print abs(det(primal_basis))
    print supposed_determinant
    # rows are basis vectors here!! That's why .T
    return primal_basis.T


print test_lattice_construction(3 ^ 5, 2, 10, 6)
