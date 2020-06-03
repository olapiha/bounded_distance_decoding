#!/usr/bin/env sage

# We assume that input data satisfies the following constrains:
# We are given a finite field F and polynomials f and g from F[x] such that deg f > deg g
# and there exist another pair of pynomials n and d that g = n/d mod f and gcd(n,d) = 1
# We also assume that deg n + deg d < 1/2 deg f
# If these conditions are satisfied there exist a unique row j of EEA(Extended Euclidean Algorithm)
# such that s_j * f + t_j * g = r_j and r_j = n, t_j = d(up to a constant)
# and this row corresponds to the unique highest degree quotient q_j with deg q_j > 1/2 deg f

# The function returns the pair (n, d) or throws
# an exception if the constrains above are not satisfied
def rational_function_reconstruction(g, f, Fx):
    assert f.degree() > g.degree()
    if g == 0:
        assert f.degree() > 3 #do I actually need this? Maybe just for coherence
        return (0,1)
    (r0, r1) = (f, g)
    (t0, t1) = (0, 1)
    q = Fx(1)
    while q.degree() <= f.degree()/2:
        (q, r) = r0.quo_rem(r1)
        (r0, r1) = (r1, r) # r = r0 - q*r1
        (t0, t1) = (t1, t0 - q*t1)
    # r0 and t0 is my answer
    assert (r0.degree() + t0.degree() < f.degree()/2)
    assert r0.gcd(t0)
    return (r0, t0)

Fx.<x> = PolynomialRing(GF(79))
g = 51*x^29 + 37*x^27 + 33*x^26 + 72*x^25 + 59*x^24 + 51*x^23 + 15*x^22 + 21*x^21 + 27*x^20 + 6*x^19 + 65*x^18 + 40*x^17 + 63*x^16 + 60*x^15 + 40*x^14 + 42*x^13 + 31*x^12 + 35*x^11 + 25*x^10 + 58*x^9 + 54*x^8 + 3*x^7 + 69*x^6 + 30*x^5 + 67*x^4 + 5*x^3 + 35*x^2 + 46*x + 67
f = x^30 + 44*x^29 + 30*x^28 + 27*x^27 + 33*x^26 + 66*x^25 + 9*x^24 + 54*x^23 + 44*x^22 + 57*x^21 + 32*x^20 + 67*x^19 + 9*x^18 + 54*x^17 + 13*x^16 + 21*x^15 + 10*x^14 + 62*x^13 + 26*x^12 + 68*x^11 + 7*x^10 + 61*x^9 + 72*x^8 + 11*x^7 + 68*x^6 + 58*x^5 + 74*x^4 + 11*x^3 + 8*x^2 + 36*x + 48

# expected result is an integer multiple of this pair:
# n =  3*x^2 + 68*x + 48
# d = x^10 + 52*x^9 + 4*x^8 + 19*x^7 + 16*x^6 + 55*x^5 + 31*x^4 + 41*x^3 + 23*x^2 + 65*x + 37

print rational_function_reconstruction(g, f, Fx)
