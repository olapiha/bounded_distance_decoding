#!/usr/bin/env sage
import numpy

n = 1000
k = floor(n/(2*log(n)))
print("k =", float(k))
w = 50
print("Ternary case of ISD")
min_cost =  binomial(n, w) / binomial(k, w)
argmin = 0
for i in range(1, w+1):
    cost = 2^i* binomial(n, w) / binomial(k, w-i)
    print(i, float(cost))
    if cost < min_cost:
        argmin = i
        min_cost = cost
print("w_2 = ", argmin, float(min_cost))

print("Ternary case of ISD+MiM")
min_cost =  2^w * binomial(n, w) +  binomial(n, w) / binomial(k, w)
argmin = 0
for i in range(1, w+1):
    cost = ( 2^(w-floor(i/2))* binomial(n, w) / binomial((n-k)/2, floor(i/2)) ) + 2^(floor(i/2)) * binomial(n, w) / (binomial((n-k)/2, floor(i/2)) * binomial(k, w-i))
    print(i, float(cost))
    if cost < min_cost:
        argmin = i
        min_cost = cost
print("w_2 =", argmin, float(min_cost))
