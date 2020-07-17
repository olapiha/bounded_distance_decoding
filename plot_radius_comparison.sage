#!/usr/bin/env sage
import numpy



def upperbound(n):
    return (n,n/2)

def algorithm_basic(n):
    m = 3^n
    B = (n+1)*ln(n+1)
    return (n, ln(m/2)/(4*euler_phi(m)^(1/n)*ln(B)))

def algorithm_general(n,k):
    m = 1
    for i in range(k):
        P = Primes()
        q = P.unrank(i+1)
        q = q**round((n*ln(3))/(ln(q)*k))
        m = m * q
    B = (n+1)*ln(n+1)
    return (n, ln(m/2)/(4*euler_phi(m)^(1/n)*ln(B)))

def polynomials(n):
    return (n, (n)/(ln(n)*(n^2-1)^(1/(2*ln(n)))))

n = 180
s = 80
k = 20

p1 = list_plot([upperbound(i) for i in range(s,n)], color='red')
p3 = list_plot([algorithm_basic(i) for i in range(s,n)], color='green', frame=True, axes_labels=['$n$','$\overline{r}_{1}$'],axes=False)
p4 = list_plot([algorithm_general(i,k) for i in range(s,n)], color='orange')
p2 = list_plot([polynomials(i) for i in range(s,n)], color='black')
p = p3+p4#+p1+p2
show(p)
