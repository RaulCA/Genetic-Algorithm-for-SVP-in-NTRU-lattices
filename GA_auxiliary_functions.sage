### GA_auxiliary_functions.sage

from numpy import random
from math import floor, ceil, sqrt

R.<x> = PolynomialRing(ZZ)

def sum_mod3_centered(x, y):
    # x, y in [-1,0,1]
    return (x + 1 + y) % 3 - 1


def round_vector(myvec):
    # tested okay
    return vector([round(v) for v in myvec])

def crossover(v1, v2):
    #v1 and v2 sage vectors
    l = len(v1)
    assert l == len(v2)
    z = zip(v1, v2)
    output = vector( (0 for _ in range(l)) )
    i = 0
    for a, b in z:
        output[i] = random.choice((a, b))
        i += 1
    return output

def mutation(v):
    # v sage vector
    l = len(v)
    for i in range(l):
        r = random.random() # r <- random uniform(0,1)
        if r  < 1/l:
            # plus 1
            v[i] = ( (v[i] + 2) % 3) - 1
        elif r < 2/l:
            # minus 1
            v[i] = (v[i] % 3) - 1
            
    return v

def cost(v, B_inv):
    #n = len(v)
    #B = B_inv**-1
    #B_norms = vector( (B[i].norm() for i in range(n)) )
    d = v * B_inv - round_vector(v * B_inv)
    
    #D = vector( (abs(D[i]) * B[i].norm()  for i in range(n)) )
    #d = sum(D)
    return d*d # previously d*d

def cost_x(x):
    e = x - round_vector(x)
    return vector((abs(w) for w in e))

def one_time_local_search(v, B_inv):
    cur_cost = cost(v, B_inv)
    n = len(v)
    v2 = copy(v)
    output = copy(v)
    
    for i in range(n):
        v2[i] = ( (v[i] + 2) % 3) - 1
        v2cost = cost(v2, B_inv)
        
        if v2cost < cur_cost:
            cur_cost = v2cost
            output = copy(v2)
        
        v2[i] = (v[i] % 3) - 1        
        v2cost = cost(v2, B_inv)
        
        if v2cost < cur_cost:
            cur_cost = v2cost
            output = copy(v2)
        
        v2[i] = v[i]
    return output

def localsearch(v, B_inv):
    # get y not zero, keep y not zero
    while True:
        new_v = one_time_local_search(v, B_inv)
        if v == new_v:
            break
        v = new_v
    return new_v
