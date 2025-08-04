### NTRU_lattices.sage

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
def keypair():
    count = 0
    for count in range(100):
        try:
            f = randompolyd1d2(df,df-1)
            f3 = invertmodprime(f,p)
            fq = invertmodpowerof2(f,q)
            break
        except:
            pass
        assert(count < 1000)
    g = randompolyd1d2(dg,dg)
    publickey = balancedmod(p * convolution(fq,g),q)
    secretkey = f,f3
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

