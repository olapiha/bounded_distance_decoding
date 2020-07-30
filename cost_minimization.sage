#!/usr/bin/env sage
import numpy

n = 400
# first dimention of the public key matrix: n-t in their notation
d = 200
k = n - d
# weight of the error(?): d-1 in their notation
w = d-1
q = next_prime(400)



n=600
d=1
k=floor((n)/(ln(n)))
w=floor((d*k)/((q**d - 1)**(k/n)))
q=next_prime(n)

print(len(bin(binomial(n-d, ((n-d)*(d-1)/n).round()))))

print("n = {} \\\\".format(float(n)))
print("k = {} \\\\".format(float(k)))
print("w = {} \\\\".format(float(w)))
print("q = {} \\\\".format(float(q)))

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
    return (argmin[0], argmin[1], float(min_cost))

print("Binary Bruteforce \\\\")

cost = binomial(n,w)*(w+1)*k
print("cost = {:e} \\\\".format(float(cost)))

print("Ternary Bruteforce \\\\")

cost = binomial(n,w)*(w+1)*k*2^w
print("cost = {:e} \\\\".format(float(cost)))

print("Binary case of ISD \\\\")

isd_cost_bin = lambda w2: (k*binomial(n, w))*(n^2 +(n-1)*(k-1)
                        + w2*binomial((n-k),w2))/(binomial(k, w-w2)*binomial((n-k), w2))

res = cost_minimization(isd_cost_bin)
print ("$w_2$, cost = {}, {:e} \\\\".format(res[0], res[1]))

print("Ternary case of ISD \\\\")

isd_cost_ter = lambda w2: (k*binomial(n, w))*(n^2 +(n-1)*(k-1)
                        + w2*binomial((n-k),w2)*2^w2)/(binomial(k, w-w2)*binomial((n-k), w2))

res = cost_minimization(isd_cost_ter)
print ("$w_2$, cost = {}, {:e} \\\\".format(res[0], res[1]))

print("Binary case of MitM \\\\")

print("cost = {:e} \\\\".format(float((k*binomial(n,w)*(n^2+w*binomial((n/2).round(),(w/2).round())))/(binomial((n/2).round(),(w/2).round())^2))))

print("Ternary case of MitM \\\\")

print("cost = {:e} \\\\".format(float((k*binomial(n,w)*2^(w/2)*(n^2+w*binomial((n/2).round(),(w/2).round())))/(binomial((n/2).round(),(w/2).round())^2))))

print("Binary case of ISD+MitM \\\\")

cost_bin = lambda w2, p: (k*binomial(n,w)*(n^2 + (n-1)*(k-1) + w2*binomial(((n-k)/2).round(),(w2/2).round())
                            + ((binomial(n,w2)*p^n)/(q^n))))/(binomial(((n-k)/2).round(),(w2/2).round())^2
                            *binomial(k,w-w2)*( 1 - (binomial(k,w-w2) * (w-w2) * binomial(((n-k)/2).round(),(w2/2).round()))/(q^(k)*p)))

res = cost_minimization2(cost_bin)
print ("$w_2$, p, cost = {}, {}, {:b} \\\\".format(res[0], res[1], floor(res[2])))

print("Ternary case of ISD+MitM \\\\")

cost_ter = lambda w2, p: (k*binomial(n,w)*(n^2 + (n-1)*(k-1) + w2*2^(w2/2)*binomial(((n-k)/2).round(),(w2/2).round())
                            + ((binomial(n,w2)*2^(w2/2)*p^n)/(q^n))))/(binomial(((n-k)/2).round(),(w2/2).round())^2
                            *binomial(k,w-w2)*( 1 - (binomial(k,w-w2) * (w-w2) * 2^(w-w2) * binomial(((n-k)/2).round(),(w2/2).round()))/(q^(k)*p)))

res = cost_minimization2(cost_ter)
print ("$w_2$, p, cost = {}, {}, {:b} \\\\".format(res[0], res[1], floor(res[2])))
