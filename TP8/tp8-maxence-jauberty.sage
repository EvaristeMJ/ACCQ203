print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP8 : PRIMALITE DES ENTIERS                                                 #
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
# TEST DE RABIN-MILLER 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n = 561

# Code pour l'EXERCICE

def get_vm(n):
    """
    Entrée : n
    Sortie : (t,m) tel que n = m*(2^t) où m est impair.
    """
    temp = n.p_primary_part(2) #on récupère la partie 2^v
    m = n//temp  #on récupère la partie première avec 2
    v = len(temp.divisors())-1 #v = (nb de diviseurs de 2^v) - 1
    assert(m*2^v == n)
    return v,m


def testRM(n):
    a = randint(1,n-1)
    v,m = get_vm(n-1)
    g = gcd(a,n)
    if g>1:
        return False
    b = (a%n)^m
    if b == 1:
        return True
    for i in range(1,v+1):
        if (b^2 % n) == 1:
            g = gcd(b+1,n)
            if g==1 or g==n:
                return True
            else:
                return False
        b = b^2 % n
    return False

# Témoin RM :

v,m = get_vm(n-1)
est_temoin = True
a =10
if gcd(a,n) != 1:
    est_temoin = False
for d in range(v):
    if ZZ(a^(m*(2^d))).mod(n) == 1:
        est_temoin = False
        break

# # Affichage des resultats

print("Test de la primalite de n=",n,"avec implementation de Rabin-Miller")
print(testRM(n))
print("{} est un témoin".format(a))
print("Résultat du test de témoin : ",est_temoin)

print("""\
# ****************************************************************************
#  PERFORMANCES DE RABIN-MILLER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

nmin=10
nmax=500
nbtests = 25

# Code pour l'EXERCICE

rep2 = "1,4,9,15,21 sont assez déclarés faussement premiers. C'est le cas de beaucoup de premiers nombres. Si l'on regarde pour des valeurs plus grandes, autour de 100-140, on remarque que les nombres qui posent problèmes sont de la forme pq avec p,q deux premiers. Assez largement, il semblerait que plus il y a de nombres premiers distincts dans la décomposition, moins il y a de chances d'erreurs. Par exemple, 133 = 7x19 121 = 11x11 125, 119 etc... sont des nombres qui ont tendance à être faussement déclarés premiers."
rep3 = "On aura au moins une probabilité de 1/4 d'obtenir aléatoirement a témoin. La probabilité de ne pas avoir de témoin en k essais sera de 4^(-k) donc il faudra au moins 25 essais pour avoir une probabilité d'au moins 2^(-50)"
# # Affichage des resultats

data = [sum( [int(ZZ(n).is_prime() == testRM(n)) for i in range(nbtests)])/nbtests for n in range(nmin,nmax)] # on affiche les cas où testRM a raison.
bar_chart(data)
print(rep2)
print(rep3)
list_plot( [timeit( 'testRM(n)', number=20, repeat=3, seconds=true) for n in range(1001,1001+100000,100) ])

reset()
from random import randint
print("""\
# ****************************************************************************
# TEST DE SOLOVAY-STRASSEN 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

n = 561

# Code pour l'EXERCICE

def testSS(n):
    if n%2 == 0:
        return False
    a = randint(2,n)
    if gcd(a,n) !=1:
        return False
    g = a^((n-1)/2)%n
    if (jacobi_symbol(a,n)%n) == g:
        return True
    
    return False

nmin=10
nmax=500
nbtests = 100

rep3 = "Il semble y avoir à nouveau des problèmes avec les nombres composés qui s'écrivent comme le produit de premiers. "
rep4 = "Les candidats pour les témoins de Solovay sont dans les inversibles de Zn donc il y en a phi(n), sur lesquels on peut réaliser un tirage aléatoire uniforme. On aura alors une probabilité d'au moins 1/2 de tomber sur un témoin ainsi. La chance de ne pas tomber sur un témoin en k essais est donc d'au moins (1/2)^k. Il faudra donc une cinquantaine d'essais pour avoir une probilité de réussite de 2^(-50). On ajoute alors une recherche d'au moins 50 potentiels témoins dans l'algorithme"


# # Affichage des resultats

print("Test de la primalite de n=",n,"avec implementation de Solovay-Strassen")
print(testSS(n))
print(rep3)
print(rep4)

bar_chart( [1/nbtests * sum( [testSS(n) for i in range(nbtests)]) for n in range(nmin,nmax)])


reset()
print("""\
# ****************************************************************************
# COMPARAISON ENTRE LES TESTS DE R-M ET S-S 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

nmax=150

# Code pour l'EXERCICE
def get_vm(n):
    """
    Entrée : n
    Sortie : (v,m) tel que n = m*(2^t) où m est impair.
    """
    temp = n.p_primary_part(2) #on récupère la partie 2^t
    m = n//temp  #on récupère la partie première avec 2
    v = len(temp.divisors())-1 #t = (nb de diviseurs de 2^t) - 1
    assert(m*2^v == n)
    return v,m

def est_temoin_RM(n,a):
    if gcd(a,n) != 1:
        return False
    v,m = get_vm(n-1)
    est_temoin = True
    for d in range(v):
        if ZZ(a^(m*(2^d))).mod(n) == 1:
            est_temoin = False
            break
    return est_temoin

def est_temoin_SS(n,a):
    if gcd(a,n)!=1:
        return False
    jac = jacobi_symbol(a,n)
    temp = a^((n-1)/2)
    return (jac-temp)%n != 0

Temoins = []
for n in range(3,150,2):
    if ZZ(n).is_prime():
        continue
    for a in range(1,n):
        if est_temoin_RM(n,a): #On vérifie si a est témoin RM
            if not est_temoin_SS(n,a): # puis on vérifie qu'il n'est pas témoin SS
                Temoins.append((n,a))

# # Affichage des resultats

print("Liste d'entiers composés et de temoins exclusifs de Rabin-Miller")
print(Temoins)



reset()
print("""\
# ****************************************************************************
# TEST DE LUCAS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

def get_tm(n):
    """
    Entrée : n
    Sortie : (t,m) tel que n = m*(2^t) où m est impair.
    """
    temp = n.p_primary_part(2) #on récupère la partie 2^t
    m = n//temp  #on récupère la partie première avec 2
    t = len(temp.divisors())-1 #t = (nb de diviseurs de 2^t) - 1
    assert(m*2^t == n)
    return t,m


def testL(n,p,q):
    if n == 2:
        return True
    
    d = p^2-4*q
    g = gcd(n,2*q*d)
    if 1<g<n:
        return False
    if g==n:
        return None
    
    t,m = get_tm(n- jacobi_symbol(d,n))
    
    # On construit la suite u
    U = [0,1]
    for i in range(2,m+1) :
        U.append(p*U[i-1]-q*U[i-2])
    
    # On construit la suite v
    V = [2,p]
    for i in range(2,n+1) :
        V.append(p*V[i-1]-q*V[i-2])
        
    g = gcd(n,U[m])
    
    if  1<g<n:
        return False
    if g==n:
        return True
    for s in range(t):
        g = gcd(n,V[m*(2^s)])
        if  1<g<n:
            return False
        if g==n:
            return True
    return False 
# # Affichage des resultats

print("Test de Lucas sur 10 essais")
for _ in range(10):
    n =  ZZ.random_element(2,250)
    print("Pour n = {}, le test de Lucas est juste".format(n),n.is_prime()==testL(n,3,14))




reset()
print("""\
# ****************************************************************************
# TEST DE BAILLIE, POMERANCE, SELFRIDGE ET WAGSTAFF
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

nmax=1000

# Code pour l'EXERCICE

def get_tm(n):
    """
    Entrée : n
    Sortie : (t,m) tel que n = m*(2^t) où m est impair.
    """
    temp = n.p_primary_part(2) #on récupère la partie 2^t
    m = n//temp  #on récupère la partie première avec 2
    t = len(temp.divisors())-1 #t = (nb de diviseurs de 2^t) - 1
    assert(m*2^t == n)
    return t,m

def testL(n,p,q):
    if n == 2:
        return True
    
    d = p^2-4*q
    g = gcd(n,2*q*d)
    if 1<g<n:
        return False
    if g==n:
        return None
    
    t,m = get_tm(n- jacobi_symbol(d,n))
    
    # On construit la suite u
    U = [0,1]
    for i in range(2,m+1) :
        U.append(p*U[i-1]-q*U[i-2])
    
    # On construit la suite v
    V = [2,p]
    for i in range(2,n+1) :
        V.append(p*V[i-1]-q*V[i-2])
        
    g = gcd(n,U[m])
    
    if  1<g<n:
        return False
    if g==n:
        return True
    for s in range(t):
        g = gcd(n,V[m*(2^s)])
        if  1<g<n:
            return False
        if g==n:
            return True
    return False 

def testRM(n,a):
    v,m = get_tm(n-1)
    g = gcd(a,n)
    if g>1:
        return False
    b = (a%n)^m
    if b == 1:
        return True
    for i in range(1,v+1):
        if (b^2 % n) == 1:
            g = gcd(b+1,n)
            if g==1 or g==n:
                return True
            else:
                return False
        b = b^2 % n
    return False
def testBPSW(n):
    if n == 2:
        return True
    if not testRM(n,2):
        return False
    d = 5
    k=0
    while jacobi_symbol(d,n) != -1:
        k= k+1
        d = ((-1)^k)*(2*k+5)

    
    p = 1
    q = (1-d)/4
    if not testL(n,p,q):
        return False
    
    return True
    

# # Affichage des resultats

print(all([ZZ(n).is_prime()==testBPSW(n) for n in range(2,nmax+1)]))

