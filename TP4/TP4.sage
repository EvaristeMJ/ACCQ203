print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP4 : FACTORISATION DE POLYNOMES UNIVARIEES SUR CORPS FINIS                 #
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
# FACTORISATION DES PUISSANCES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F1849.<omega> = FiniteField(43^2,modulus=x^2+1)
Pol1849.<x> = PolynomialRing(F1849)
f=x^172+(3-2*omega)*x^129-5*omega*x^86+(2 + 4*omega)*x^43-1-omega 

F9.<alpha> = FiniteField(9)
Pol9.<y> = PolynomialRing(F9)
g = y^30-y^15+alpha*y^3+1

# Code pour l'EXERCICE

def racine_p_polynome(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    u=Pol(0)
    l_pol = f.list()
    for i in range(len(l_pol)):
        if l_pol[i] != 0:
            d = i//p
            c = l_pol[i].pth_root()
            u = u + c*x^d
    assert(u^p==f)
    return u

test = true
n_test = 100

for i in range(n_test):
    h = Pol1849.random_element()
    if h!= racine_p_polynome(h^43):
        test = false
        break


# # Affichage des resultats

print( "\n$ Question 3")
print( "La racine de",f,"est",racine_p_polynome(f))
print( "\n$ Question 4")
print( "Test sur 100 exemples : ",test)



reset()
print("""\
# ****************************************************************************
# FACTORISATION SANS FACTEURS CARRES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F7 = FiniteField(7)
Pol7.<x> = PolynomialRing(F7)
f = x^10 +6*x^9 +3*x^7 +3*x^3 +4*x^2 +2

# Code pour l'EXERCICE

def myFsFC(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [(f,1)]
    if f.degree()<=0:
        return []
    if f.derivative() !=0:
        i = 1
        retour = []
        t =  xgcd(f,f.derivative())[0]
        u = f//t
        while u != 1:
            y = xgcd(t,u)[0]
            if i%p!=0 and u//y != 1:
                retour.append((u//y,i))
            i = i +1
            u = y
            t = t//y
        if t != 1:
            K = myFsFC(t.nth_root(p))
            for i in range(len(K)):
                K[i] = (K[i][0],K[i][1]*p)
            retour = retour + K
    else:
        K = myFsFC(f.nth_root(p))
        for i in range(len(K)):
            K[i] = (K[i][0],K[i][1]*p)
        retour = retour + K
            
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour


test = true
n_test = 1000

def avec_facteurs_carres(h):
    if h.derivative() == 0:
        return true
    if xgcd(h.derivative(),h)[0] != 1:
        return true
    return false

for i in range(n_test):
    h = Pol7.random_element(degree = (5,10))
    L = myFsFC(h)
    for (g,e) in L:
        if g.degree() <= 0:
            continue
        if avec_facteurs_carres(g):
            test = false
# # Affichage des resultats

print( "\n$ Question 2")
print( "La factorisation de",f,"est",myFsFC(f))
print( "\n$ Question 4")
print( "Test sur 100 exemples : ",test)


reset()
import random
print("""\
# ****************************************************************************
# FACTORISATION ETAGEE EN DEGRES DISTINCTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F5 = FiniteField(5)
Pol.<x>=PolynomialRing(F5)
f = x^10-2*x^9+x^8+x^7-x^6-2*x^5+2*x^4+2*x^3-x

# Code pour l'EXERCICE


def myFEDD(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    h = x
    g = 1
    fi = f
    while fi!=1:
        h = h^q
        h = h.mod(f)
        g = xgcd(h-x,fi)[0]
        fi = fi//g
        retour.append(g)
    assert(prod(retour) == f)
    return retour
    
L1 = [f for f in Pol.polynomials(of_degree=1) if f.is_irreducible() and f.is_squarefree()
 and f.leading_coefficient()==1]
L2 = [f for f in Pol.polynomials(of_degree=2) if f.is_irreducible() and f.is_squarefree()
 and f.leading_coefficient()==1]
L3 = [f for f in Pol.polynomials(of_degree=3) if f.is_irreducible() and f.is_squarefree()
 and f.leading_coefficient()==1]
test = true
n_test = 100
for i in range(n_test):
    h1 = L1[random.randint(0,len(L1)-1)]
    h2 = L2[random.randint(0,len(L2)-1)]
    h3 = L3[random.randint(0,len(L3)-1)]
    g=h1*h2*h3
    fedd = myFEDD(g)
    if not (h1 in fedd and h2 in fedd and h3 in fedd):
        test = false
        break



# # Affichage des resultats

print( "\n$ Question 1")
print( "La factorisation de",f,"est",myFEDD(f))
print("Résultat du test de la fonction myFEDD sur 100 exemples")
print(test)

reset()
import random
from collections import Counter
print("""\
# ****************************************************************************
# CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=3
d=4
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 


# Code pour l'EXERCICE

def Tr2(m,f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    tr = Pol(0)
    for i in range(m):
        tr = tr+f^(2^i)
    return tr
def myCZ(f,d):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    if f.degree() == d:
        return retour
    u = Pol.random_element(degree = (0,2*d-1))
    b = xgcd(f,u^(q^d))[0]
    while b.degree()<=0 or b.degree()>= f.degree():
            u = Pol.random_element(degree = (0,2*d-1))
            b = xgcd(f,u^(q^d))[0]
    retour1 = myCZ(b,d)
    retour2 = myCZ(f//b,d)
    retour = retour1 + retour2
    assert(prod(retour) == f)
    return retour
    
def myCZ2(f,d):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    k = int(q.log(p))
    retour = [f]
    if f.degree() == d:
        return retour
    else:
        u = Pol.random_element(degree = (0,2* d-1))
        tr = Tr2(k,u)
        b = xgcd(f,tr)[0]
        while b.degree()<=0 or b.degree>= f.degree():
            u = Pol.random_element(degree = (0,2*d-1))
            tr = Tr2(k,u)
            b = xgcd(f,tr)[0]
        retour1 = myCZ2(b,d)
        retour2 = myCZ2(f//b,d)
        retour = retour1 + retour2
    assert(prod(retour) == f)
    return retour    

test  = true
n = 100
L =[f for f in Polq.polynomials(of_degree=d) if f.is_irreducible()
 and f.leading_coefficient()==1]
n_facteurs = 10
for i in range(n):
    m = random.randint(2,n_facteurs)
    K = random.sample(L,m)
    f = prod(K)
    facteurs = myCZ(f,d)
    if Counter(K) != Counter(facteurs): # vérifie si les facteurs ne sont pas les mêmes
        test = false
        break

F2=FiniteField(2)
Pol2.<y> = PolynomialRing(F2)

test2 = true
L2 = [f for f in Pol2.polynomials(of_degree=d) if f.is_irreducible()
 and f.leading_coefficient()==1]

for i in range(n):
    m = random.randint(2,len(L2))
    K2 = random.sample(L2,m)
    f = prod(K2)
    facteurs = myCZ(f,d)
    if Counter(K2) != Counter(facteurs):
        test2 =false
        break
# # Affichage des resultats
print("Résultat du test sur 100 polynômes aléatoires de myCZ :")
print(test)
print("Résultat du test sur 100 polynômes aléatoires de myCZ2 :")
print(test2)


reset()
print("""\
# ****************************************************************************
# FACTORISATION COMPLETE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]
    
f = L1[0]*L1[1]^3*L1[2]^4
g= L2[0]*L2[1]^4*L2[2]^4
h= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4


# Code pour l'EXERCICE

def myFsFC(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [(f,1)]
    if f.degree()<=0:
        return []
    if f.derivative() !=0:
        i = 1
        retour = []
        t =  xgcd(f,f.derivative())[0]
        u = f//t
        while u != 1:
            y = xgcd(t,u)[0]
            if i%p!=0 and u//y != 1:
                retour.append((u//y,i))
            i = i +1
            u = y
            t = t//y
        if t != 1:
            K = myFsFC(t.nth_root(p))
            for i in range(len(K)):
                K[i] = (K[i][0],K[i][1]*p)
            retour = retour + K
    else:
        K = myFsFC(f.nth_root(p))
        for i in range(len(K)):
            K[i] = (K[i][0],K[i][1]*p)
        retour = retour + K
            
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour


def myFEDD(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    h = x
    g = 1
    fi = f
    while fi!=1:
        h = h^q
        h = h.mod(f)
        g = xgcd(h-x,fi)[0]
        fi = fi//g
        retour.append(g)
    assert(prod(retour) == f)
    return retour

def Tr2(m,f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    tr = Pol(0)
    for i in range(m):
        tr = tr+f^(2^i)
    return tr
def myCZ(f,d):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    if f.degree() == d:
        return retour
    u = Pol.random_element(degree = (0,2*d-1))
    b = xgcd(f,u^(q^d))[0]
    while b.degree()<=0 or b.degree()>= f.degree():
            u = Pol.random_element(degree = (0,2*d-1))
            b = xgcd(f,u^(q^d))[0]
    retour1 = myCZ(b,d)
    retour2 = myCZ(f//b,d)
    retour = retour1 + retour2
    assert(prod(retour) == f)
    return retour

def myFactorisation(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    for (g,e) in myFsFC(f):
        G = myFEDD(g)
        for i in range(0,len(G)):
            if G[i] == 1:
                continue
            for fi in myCZ(G[i],i+1):
                retour.append((fi,e))
    
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour

facteurs_f = Factorization(myFactorisation(f))
facteurs_f.simplify()
facteurs_g = Factorization(myFactorisation(g))
facteurs_g.simplify()
facteurs_h = Factorization(myFactorisation(h))
facteurs_h.simplify()

test = facteurs_f == factor(f) and facteurs_h == factor(h) and facteurs_g == factor(g)

# # Affichage des résultats


print("Le polynôme f = ",f," se factorise en ")
print(facteurs_f)
print()
print("Le polynôme g = ",g," se factorise en ")
print(facteurs_g)
print()
print("Le polynôme h = ",h," se factorise en ")
print(facteurs_h)
print()
print("Résultat du test sur f,g,h ")
print(test)


reset()
print("""\
# ****************************************************************************
# RACINES D'UN POLYNOME
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]
    
f = L1[0]*L1[1]^3*L1[2]^4
f *= L2[0]*L2[1]^4*L2[2]^4
f *= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4


# Code pour l'EXERCICE

def myCZ(f,d):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    if f.degree() == d:
        return retour
    u = Pol.random_element(degree = (0,2*d-1))
    b = xgcd(f,u^(q^d))[0]
    while b.degree()<=0 or b.degree()>= f.degree():
            u = Pol.random_element(degree = (0,2*d-1))
            b = xgcd(f,u^(q^d))[0]
    retour1 = myCZ(b,d)
    retour2 = myCZ(f//b,d)
    retour = retour1 + retour2
    assert(prod(retour) == f)
    return retour



def myRacine(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    r = x^q-x
    r = r.mod(f)
    g = xgcd(r,f)[0]
    facteurs = myCZ(g,1)
    retour =[]
    for gi in facteurs:
        retour.append(gi.coefficients()[0])
    assert(f(z)==0 for z in retour)
    return retour

# # Affichage des resultats

print( "\n$ Question 1")
print( "Les racines de ",f,"sont",myRacine(f))

reset()
print("""\
# ****************************************************************************
# ETUDE DE CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice
d = 3
q = 7
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 

L = [f for f in Polq.polynomials(of_degree=d) if f.is_irreducible()
and f.leading_coefficient()==1]



f1 = L[1]
f2 = L[7]
# Code pour l'EXERCICE

# g est un carré modulo f_i non nul si et seulement si g^{(q^d-1)/2} = 1 mod f_i

def est_carre_mod(g,f):
    if g == 0:
        return true
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    d=f.degree()
    h = g^((q^d-1)//2)
    return h.mod(f) == 1

L1 = []
L2 = []
K1 = []
K2 = []

for g in Polq.polynomials(max_degree = d):
    if est_carre_mod(g,f1):
        g = g.mod(f1)
        if g not in L1:
            L1.append(g)
    else:
        g = g.mod(f1)
        if g not in K1:
            K1.append(g)
    if est_carre_mod(g,f2):
        g = g.mod(f2)
        if g not in L2:
            L2.append(g)
    else:
        g = g.mod(f2)
        if g not in K2:
            K2.append(g)

# # Affichage des resultats
