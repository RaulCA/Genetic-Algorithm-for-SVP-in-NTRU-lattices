### Adaptation from LatticeHacks to generate NTRU Lattices

from numpy import random
from math import floor, ceil, sqrt


global_poly_f = None
global_poly_g = None


Zx.<x> = ZZ[]
p = 3
# Cyclic convolution. This is the multiplication operation used in NTRU
def convolution(f,g):
      return (f * g) % (x^n-1)

def balancedmod(f,q):
  g = list(((f[i] + q//2) % q) - q//2 for i in range(n))
  return Zx(g)


def randomdpoly(d):
  assert d <= n
  result = n*[0]
  for j in range(d):
    while True:
      r = randrange(n)
      if not result[r]: break
    result[r] = 1-2*randrange(2)
  return Zx(result)

def randompolyd1d2(d1,d2):
  assert d1 + d2 <= n
  result = n*[0]
  pos = sample(range(n),d1+d2)
  for i in pos[0:d1]:
    result[i] = 1
  for i in pos[d1:(d1+d2)]:
    result[i] = -1
  return Zx(result)

def polynorm(a):
    return(max(a.coefficients())-min(a.coefficients()))

def polynorm2(a):
    return(max(map(abs,a.coefficients())))


def invertmodprime(f,p):
    T = Zx.change_ring(Integers(p)).quotient(x^n-1)
    return Zx(lift(1 / T(f)))

# Division modulo powers of 2
def invertmodpowerof2(f,q):
    assert q.is_power_of(2)
    g = invertmodprime(f,2)
    while True:
        r = balancedmod(convolution(g,f),q)
        if r == 1: return g
        g = balancedmod(convolution(g,2 - r),q)
        
# NTRU key generation
# remember is fixed, p = 3.
def keypair():
    count = 0
    #for count in range(10000):
    while True:
        try:
            f = randompolyd1d2(df,df-1)
            f3 = invertmodprime(f,p)
            fq = invertmodpowerof2(f,q)
            break
        except:
            pass
        #assert(count < 1000)
    
    g = randompolyd1d2(dg,dg)
    publickey = balancedmod(p * convolution(fq,g),q)
    secretkey = f,f3
    global global_poly_f
    global global_poly_g
    global_poly_f = f
    global_poly_g = g
    print('(g,f) = ', (g,f))
    return publickey,secretkey

def ntru_lattice(n_local = 10,df_local = 3, dg_local = 3, dphi_local = 3, q_local = 256):
    global dg,df,dphi,n,q
    dg = dg_local
    df = df_local
    dphi = dphi_local
    n = n_local
    q = q_local
    publickey,secretkey = keypair()
    recip3 =  lift(1/Integers(q_local)(3))
    publickeyover3 = (recip3 * publickey) % q
    M = matrix(2 * n)
    for i in range(n):
        M[i,i] = q
    for i in range(n):
        M[i+n,i+n] = 1
        c = convolution(x^i,publickeyover3)
        for j in range(n):
            M[i+n,j] = c[j]
    return M
   

R.<x> = PolynomialRing(ZZ)

#### Auxiliary functions to Genetic Algorithm for NTRU Lattices
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

#def cost_x(x):
#    e = x - round_vector(x)
#    return vector((abs(w) for w in e))

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
    
    
 
B = ntru_lattice(n_local = 50, df_local = 13, dg_local = 13, dphi_local = 12, q_local = 256)
B_inv = B**-1; print(B_inv)


# complete the vectors with 0's if necessary (the case polynomial has degree less than n-1)
poly_g = vector(global_poly_g).list()
poly_f = vector(global_poly_f).list()
if len(poly_g) < n:
    poly_g = poly_g + [0 for i in range(n - len(poly_g))]
if len(poly_f) < n:
    poly_f = poly_f + [0 for i in range(n - len(poly_f))]

s = vector(poly_g + poly_f); print('s = ', s); print('len(s) = ', len(s)); print('s*B_inv = ', s*B_inv)
n = B.nrows()



### Generate sample cost values with and without BKZ-reduction
sample_size = 10000
matrix_cost_neigh_no_bkz  = matrix(QQ, [[-1 for i in range(1, n)] for j in range(sample_size)])
matrix_cost_neigh_yes_bkz = matrix(QQ, [[-1 for i in range(1, n)] for j in range(sample_size)])
print( matrix_cost_neigh_yes_bkz.nrows(), matrix_cost_neigh_yes_bkz.ncols() )

### costs WITHOUT and WITH BKZ for multiple Hamming distances
for ham_dist in range(1, n):
    neigh1 = matrix([solution for i in range(sample_size)])
    # the 'sample' method returns a list object. Using sample(tuple/list , n)[0] to get the first element of the output
    for i in range(sample_size):
        rvec = sample(range(n), ham_dist)
        for r in rvec:
            neigh1[i, r] = (neigh1[i, r] + 1 + sample((-1, 1), 1)[0]) % 3 - 1
        #print('i = ', i)
        matrix_cost_neigh_no_bkz[i, ham_dist - 1]  = cost(neigh1[i], B_inv)
        matrix_cost_neigh_yes_bkz[i, ham_dist - 1] = cost(neigh1[i], B_BKZ_inv)
        


### Generate violin plots relation Hamming distance from the key to cost values    
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = list(matrix_cost_neigh_no_bkz.transpose())
data = data[0:50]

# Create sample data
np.random.seed(0)
df = pd.DataFrame({
    'Cost Value': np.concatenate(data),
    'Hamming Distance': [str(x) for x in range(1, 51) for i in range(sample_size)]
})

# Plot multiple violin plots
plt.figure(figsize=(35, 8))
sns.violinplot(x='Hamming Distance', y='Cost Value', data=df, inner='box', palette='Set2')
plt.title("Using Original Basis: Distribution of Hamming Distances X Cost Values")
plt.grid(True)
plt.show()

data = list(matrix_cost_neigh_yes_bkz.transpose())
data = data[0:50]

# Create sample data
np.random.seed(0)
df = pd.DataFrame({
    'Cost Value': np.concatenate(data),
    'Hamming Distance': [str(x) for x in range(1, 51) for i in range(sample_size)]
})

# Plot multiple violin plots
plt.figure(figsize=(35, 8))
sns.violinplot(x='Hamming Distance', y='Cost Value', data=df, inner='box', palette='Set2')
plt.title("Using BKZ-reduced Basis: Distribution of Hamming Distances X Cost Values")
plt.grid(True)
plt.show()
