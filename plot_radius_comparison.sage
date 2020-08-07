#!/usr/bin/env sage
import numpy



def upperbound(n):
    return (n,n/(e*2))

def algorithm_basic(n):
    m = 3^n
    B = (n+1)*ln(n+1)
    return (n, ln(m/2)/(4*euler_phi(m)^(1/n)*ln(B)))

def algorithm_general_optimal(n):
    m = 1
    B = (n+1)*ln(n+1)
    i = 0
    max_radius = 0
    max_m = m
    P = Primes()
    while True:
        q = P.unrank(n+i)
        m = m * q
        radius = ln(m/2)/(4*euler_phi(m)^(1/n)*ln(B))
        if radius > max_radius:
            max_radius = radius
            max_m = m
            i+=1
        else: break
    return (n, max_radius)

def polynomials(n):
    return (n, (n)/(2*ln(n)*(n^2-1)^(1/(2*ln(n)))))

n = 300
s = 10

p1 = list_plot([upperbound(i) for i in range(s,n)], color='red')
p3 = list_plot([algorithm_basic(i) for i in range(s,n)], color='green', frame=True, axes_labels=['$n$','$\overline{r}_{1}$'],axes=False)
p2 = list_plot([polynomials(i) for i in range(s,n)], color='black')
p5 = list_plot([algorithm_general_optimal(i) for i in range(s,n)], color='red')

p = p5+p3+p1+p2
show(p)
