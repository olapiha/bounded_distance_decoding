#!/usr/bin/env sage
import numpy



def upperbound(n):
    return (n,n/2)

def algorithm_basic(n):
    m = 3^n
    B = (n+1)*ln(n+1)
    return (n, ln(m/2)/(euler_phi(m)^(1/n)*ln(B)))

def algorithm_general(n,k):
    m = 1
    for i in range(k):
        P = Primes()
        q = P.unrank(i+1)
        q = q**round((n*ln(3))/(ln(q)*k))
        m = m * q
    B = (n+k)*ln(n+k)
    return (n, ln(m/2)/(euler_phi(m)^(1/n)*ln(B)))

def polynomials(n):
    return (n, (n)/(ln(n)*(n^2-1)^(1/(2*ln(n)))))

n = 100
k = 5

p1 = list_plot([upperbound(i) for i in range(10,n)], color='red')
p3 = list_plot([algorithm_basic(i) for i in range(10,n)], color='green')
p4 = list_plot([algorithm_general(i,k) for i in range(10,n)], color='orange')
p6 = list_plot([polynomials(i) for i in range(10,n)], color='black')
p = p1+p6+p3+p4#+p5+p6
show(p)
