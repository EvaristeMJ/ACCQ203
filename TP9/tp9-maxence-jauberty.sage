print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP9 : FACTORISATION DES ENTIERS                                             #
# *************************************************************************** #
# *************************************************************************** #
""")

# CONSIGNES
#
# Les seules lignes a modifier sont annoncee par "Code pour l'exercice"
# indique en commmentaire et son signalees
# Ne changez pas le nom des variables
#
# CONSEILS
#
# Ce modele vous sert a restituer votre travail. Il est deconseille d'ecrire
# une longue suite d'instruction et de debugger ensuite. Il vaut mieux tester
# le code que vous produisez ligne apres ligne, afficher les resultats et
# controler que les objets que vous definissez sont bien ceux que vous attendez.
#
# Vous devez verifier votre code en le testant, y compris par des exemples que
# vous aurez fabrique vous-meme.
#


reset()

print("""\
# ****************************************************************************
# DIVISEURS SUCCESSIFS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=2*3*3*5*5*5*7*11*11

# Code pour l'EXERCICE

def div_successives(n):
    D = []
    d = 2
    while d^2 <= n:
        while n%d == 0:
            n = n//d
            if d not in D:
                D.append(d)
        d = d+1
    if n != 1:
        D.append(n)
    return  D

# # Affichage des resultats

print("Résultat divisions successives ",div_successives(n))
for n in range(2,10):
    assert(div_successives(ZZ(n))==ZZ(n).prime_divisors())

reset()

print("""\
# ****************************************************************************
# FACTORISATION D'UN NOMBRE B-FRIABLE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n=2*3*3*5*5*5*7*11*11
P=[p for p in primes(12)]

# Code pour l'EXERCICE

def div_successives_friable(n, P):
    L = []
    for p in P:
        while n%p == 0:
            n = n//p
            if p not in L:
                L.append(p)
    if n == 1:
        return L
    
    return  "{} n'est pas friable.".format(n)

# # Affichage des resultats

print("Résultat division successible friable ",div_successives_friable(n,P))
for n in range(2,10):
    assert(div_successives_friable(ZZ(n),P)==ZZ(n).prime_divisors())

reset()
print("""\
# ****************************************************************************
# RHO DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=222763

# Code pour l'EXERCICE

def myPollardrho(n,k_b=1000):
    x0=Integers(n).random_element()
    y = x0
    def f(z):
        return z^2+1
    g = 0
    k_b = 0
    while g <= 1:
        x0 = f(x0)
        y = f(f(y))
        g = gcd(x0-y,n)
        k_b = k_b  + 1
        if k_b == 1000: # empêche le while de durer infiniment
            return 1
    if g == n:
        return 1
    return g

# # Affichage des resultats

print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardrho(n), 
      " | n est-il composé ?",not n.is_prime())

for _ in range(10):
    n=ZZ.random_element(3,100)
    print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardrho(n), 
      " | n est-il composé ?",not n.is_prime())

reset()
print("""\
# ****************************************************************************
# P-1 DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=1323269

# Code pour l'EXERCICE

def myPollardpm1(n,P):
    a = Integers(2,n).random_element()
    g = gcd(a,n)
    if g >1:
        return g
    for p in P:
        alpha = 0
        while p^alpha <= P[-1]:
            alpha = alpha + 1
        alpha = alpha -1
        
        for j in range(1,alpha + 1):
            a = Mod(a^p,n)
    
    g = gcd(a-1,n)
    if 1< g < n:
        return g
    return 1

P = [p for p in primes(20)]

# # Affichage des resultats

print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardpm1(n,P), 
      " | n est-il composé ?",not n.is_prime())

for _ in range(5):
    n=ZZ.random_element(3,100)
    print(n, 
      "| Resultat rho de Pollard : ", 
      myPollardpm1(n,P), 
      " | n est-il composé ?",not n.is_prime())

reset()
print("""\
# ****************************************************************************
# CRIBLE QUADRATIQUE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n=2886

# Code pour l'EXERCICE
def L(n):
    return exp((ln(n)*ln(ln(n)))^(1/2))

def est_friable(n,P):
    """
    entrée : liste de premiers
    sortie : True si les diviseurs premiers de n sont dans P, faux sinon
    """
    F = []
    m = n
    for i in range(len(P)):
        e = 0
        p = P[i]
        while Mod(m,p) == 0:
            m = ZZ(m)
            m = m/p
            e = e + 1
        F.append((p,e))
    return (m == 1,F)

def cribleQuadratique (n):
    B = ceil(L(n))
    P = [p for p in primes(B)]
    m = len(P)
    S = []
    c = 0
    x = ceil(sqrt(n))
    a = Mod(x^2,n)
    while c < m+1:
        b,F = est_friable(a,P)
        if b:
            S.append((x,a,F))
            c = c +1
        x = x+1
        a = Mod(x^2,n)
    E = matrix(ZZ,m,m+1)
    for i in range(m):
        for j in range(m+1):
            E[i,j] = S[j][2][i][1]
        
    E2 = matrix(FiniteField(2),E)
    v = E2.right_kernel().basis()[0]
    K = [i for i in range(len(v)) if v[i]!= 0]
    z = 1
    y = 1
    for k in K:
        x,a,_ = S[k]
        z = z*x % n
    
    for i in range(m):
        s = 0
        for k in K:
            s = s + E[i,k]
        s = s/2
        y = (y * P[i]^s) % n
        
    d = gcd(z-y,n)
    
    return  d

# # Affichage des resultats

print("Résultat crible quadratique : ",cribleQuadratique (n))