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
def parity_check_representation(q, primes):
    Zq = Zmod(q)
    multiplicative_group = Zq.unit_group()
    #Q: Should I check if it is actually cyclic?
    gen = multiplicative_group.gens_values()[0]
    logs = []
    for p in primes:
        log = compute_log(multiplicative_group, gen, p)
        logs.append(log)
    return logs

def compute_log(group, generator, element):
    pass
    #log = group.exp(element)
    #print generator**log == element
    #return log


def dual_generating_set(arg):
    pass


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


#test_pbfd(2)
#test_lll(3)

#q = 25
#Zq = Zmod(q)
#multiplicative_group = UnitGroup(Zq)
#gen = multiplicative_group.gens_values()[0]
#compute_log(multiplicative_group, gen, 7)
