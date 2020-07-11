#!/usr/bin/env sage
import numpy

n = 1000
k = floor(n/(2*log(n)))
print("k =", float(k))
w = 50
q = 50

def cost_minimization(cost_function):
    min_cost = cost_function(0)
    argmin = 0
    for i in range(1, w+1):
        cost = cost_function(i)
        #print(i, float(cost))
        if cost < min_cost:
            argmin = i
            min_cost = cost
    return (argmin, float(min_cost))

def cost_minimization2(cost_function):
    min_cost = cost_function(0,1)
    argmin = (0,1)
    for i in range(0, w):
        for j in range(1,q):
            cost = cost_function(i,j)
            #print(i, j, float(cost))
            if cost < min_cost:
                argmin = (i,j)
                min_cost = cost
    return (argmin, float(min_cost))


print("Binary case of ISD")

isd_cost_bin = lambda w2: (k*binomial(n, w))*(n^2 +(n-1)*(k-1)
                        + w2*binomial((n-k),w2))/(binomial(k, w-w2)*binomial((n-k), w2))

print ("w_2, cost = ", cost_minimization(isd_cost_bin))

print("Ternary case of ISD")

isd_cost_ter = lambda w2: (k*binomial(n, w))*(n^2 +(n-1)*(k-1)
                        + w2*binomial((n-k),w2)*2^w2)/(binomial(k, w-w2)*binomial((n-k), w2))

print ("w_2, cost = ", cost_minimization(isd_cost_ter))

print("Binary case of MitM")

print("cost =", float((k*binomial(n,w)*(n^2+w*binomial(n/2,w/2)))/(binomial(n/2,w/2)^2)))

print("Ternary case of MitM")

print("cost =", float((k*binomial(n,w)*2^(w/2)*(n^2+w*binomial(n/2,w/2)))/(binomial(n/2,w/2)^2)))

print("Binary case of ISD+MitM")

cost_bin = lambda w2, p: (k*binomial(n,w)*(n^2 + (n-1)*(k-1) + w2*binomial(((n-k)/2).round(),(w2/2).round())
                            + ((binomial(n,w2)*p^n)/(q^n))))/(binomial(((n-k)/2).round(),(w2/2).round())^2
                            *binomial(k,w-w2)*( 1 - (binomial(k,w-w2) * (w-w2) * binomial(((n-k)/2).round(),(w2/2).round()))/(q^(k)*p)))

print ("w_2, p, cost = ", cost_minimization2(cost_bin))

print("Ternary case of ISD+MitM")

cost_ter = lambda w2, p: (k*binomial(n,w)*(n^2 + (n-1)*(k-1) + w2*2^(w2/2)*binomial(((n-k)/2).round(),(w2/2).round())
                            + ((binomial(n,w2)*2^(w2/2)*p^n)/(q^n))))/(binomial(((n-k)/2).round(),(w2/2).round())^2
                            *binomial(k,w-w2)*( 1 - (binomial(k,w-w2) * (w-w2) * 2^(w-w2) * binomial(((n-k)/2).round(),(w2/2).round()))/(q^(k)*p)))

print ("w_2, p, cost = ", cost_minimization2(cost_ter))
