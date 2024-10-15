print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP5 : FACTORISATION COMPLETE DE POLYNOMES UNIVARIEES                        #
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
import random
print("""\
# ****************************************************************************
# BERLEKAMP
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F3 = FiniteField(3)
Pol3.<x> = PolynomialRing(F3)
f = x^3 - x^2 - 1

# Code pour l'EXERCICE

x3 = (x^3).mod(f)
x6 = (x^6).mod(f)


L = [f for f in Pol3.polynomials(max_degree = 3) if f.is_squarefree()
 and f.leading_coefficient()==1 and f.degree()>=1]

def myFsFC(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
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
# Pour générer la matrice de Petr Berlekamp de f
def PBmatrix(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    F = FiniteField(q)
    n = f.degree()
    Q = matrix(F,n)
    for i in range(n):
        g = (x^i)^q
        g = g.mod(f)
        for j in range(n):
            Q[i,j] = g[j]
    return Q.T

def myB(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    F = FiniteField(q)
    n = f.degree()
    
    K = PBmatrix(f)
    K = K - matrix.identity(F,n)
    
    ker = K.right_kernel()
    basis = ker.basis()
    base = []
    retour = [f]
    for vec in ker.basis():
        b = Pol3(0)
        for i in range(0,vec.length()):
            b += vec[i]*x^(i)
        base.append(b)
        
    j=0
    while len(retour) < len(base):
        j = j +1
        C = [f_tilde for f_tilde in retour if f_tilde.degree()>1]
        for f_tilde in C:
            B = []
            for alpha in F:
                a = xgcd(f_tilde,base[j]-alpha)[0]
                if a.degree() >=1:
                    B.append(a)
            retour.remove(f_tilde)
            retour = retour + B
    assert(Set(retour)== Set(g for g,_ in list(f.factor())))
    return retour

def myFactor(f):
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
            for fi in myB(G[i]):
                retour.append((fi,e))
    assert(Set(retour) == Set(list(f.factor())))
    return retour

Q=PBmatrix(f)
K = Q-matrix.identity(F3,3)

# Calcul du noyau de Q-I
ker= K.right_kernel()
b1 = ker.basis()[0]
b2 = ker.basis()[1]

n_test = 100
test = true
#Si le programme ne s'arrête pas, les conditions de sortie sont vérifiées
                                    #ce qui justifie que le code marche sur ces exemples

#Test de Berlekamp
for i in range(n_test):  
    h = L[random.randint(0,len(L)-1)]
    facteur_h = myB(h)

L2 = [f for f in Pol3.polynomials(max_degree = 5) if f.leading_coefficient()==1 and f.degree()>=2]

#Test de la fonction myFactor
for i in range(n_test):
    h = L2[random.randint(0,len(L)-1)]
    facteur_h = myFactor(h)
    
# # Affichage des resultats

print("\n$1a/ x^3 vaut",x3," et x^6 vaut",x6)
print("La matrice de Petr Berlekamp est")
print(Q)

print("\n$1b/ On a Q * b1 - b1 = ")
print(Q*b1-b1)
print("et Q * b2 - b2 = ")
print(Q*b2-b2)

print("Test de myFactor sur 100 essais :")
print(test)


reset()
print("""\
# ****************************************************************************
# RELEVEMENT DE HENSEL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

PolZZ.<x> = PolynomialRing(ZZ)
m = 5
f = x^4-1
g = x^3+2*x^2-x-2
h = x-2
d,ss,tt = xgcd(g,h)
s=PolZZ(ss/mod(d,m)); t=PolZZ(tt/mod(d,m))

# Code pour l'EXERCICE



def polynomeCentre(f,m):
    Pol=f.parent()
    x=Pol.gen()
    coeff = f.list()
    f_tilde = Pol(0)
    for i in range(len(coeff)):
        c = coeff[i].mod(m)
        if c > m//2:
            c = c-m
        f_tilde =f_tilde + c*x^i
    retour = f_tilde
    return retour

def myHensel(f,g,h,s,t,m):
    Pol=f.parent()
    x=Pol.gen()
    e = polynomeCentre(f-g*h,m^2)
    q,r = (s*e).quo_rem(h)
    q,r = polynomeCentre(q,m^2),polynomeCentre(r,m^2)
    g_etoile = polynomeCentre(g + t*e + q*g,m^2)
    h_etoile = polynomeCentre(h+r,m^2)
    b = polynomeCentre(s*g_etoile +t*h_etoile-1,m^2)
    c,d = (s*b).quo_rem(h_etoile)
    c,d = polynomeCentre(c,m^2),polynomeCentre(d,m^2)
    s_etoile = polynomeCentre(s-d,m^2)
    t_etoile = polynomeCentre(t - t*b-c*g_etoile,m^2)

    retour = g_etoile,h_etoile,s_etoile,t_etoile
    return retour

def myHenselItere(f,g,h,s,t,m,l):
    Pol=f.parent()
    x=Pol.gen()
    d,ss,tt = xgcd(g,h)
    s=Pol(ss/mod(d,m)); t=PolZZ(tt/mod(d,m))
    i = m
    ff = f
    gg=g
    hh = h
    ss=s
    tt = t
    while i < m^l:
        gg,hh,ss,tt = myHensel(ff,gg,hh,ss,tt,i)
        i = i^2
    retour = gg,hh
    return retour

reponseQ5="Soit f un polynôme de Z[X]. Si jamais f = g_1*...*g_(n+1) avec les gi premiers entre eux, donc en particulier h_n = g1*...*g_n et g_{n+1} sont premiers, on peut utiliser l'algorithme pour avoir hh_n et gg_(n+1) tel que f=hh_n * gg_(n+1). On a ensuite que hh =g1*...*g_n mod p^l, avec hh = (g_1*...*g_(n-1))*g_n mod p, on peut alors procéder par récurrence pour déterminer, gg_1,...,gg_(n+1) tel que f = gg_1*...*gg_(n+1) mod p^l"
test = false

# # Affichage des resultats

print("\n$1b/ Relèvement de ",f,"= (",g,")*(",h,")")
print(myHensel(f,g,h,s,t,m))
g_etoile,h_etoile,s_etoile,t_etoile = myHensel(f,g,h,s,t,m)
def test_myHensel(f,g,g_etoile,h,h_etoile,s,s_etoile,t,t_etoile,m):
    test = true
    test = polynomeCentre(g_etoile,m)== polynomeCentre(g,m) and polynomeCentre(h_etoile,m)== polynomeCentre(h,m) and polynomeCentre(s_etoile,m)== polynomeCentre(s,m) and polynomeCentre(t_etoile,m)== polynomeCentre(t,m)
    test = test and s_etoile.degree()<h_etoile.degree() and t_etoile.degree()<g_etoile.degree()
    test = test and polynomeCentre(f,m^2)== polynomeCentre(g_etoile*h_etoile,m^2) and polynomeCentre(s_etoile*g_etoile+t_etoile*h_etoile,m^2)== 1
    return test

test = test_myHensel(f,g,g_etoile,h,h_etoile,s,s_etoile,t,t_etoile,m)
print("Vérification des hypothèses (8) de sortie de myHensel")
print(test)
print("")
print("myHenselItere modulo 25")
print(myHenselItere(f,g,h,s,t,m,2))
print("A-t-on f-gg*hh = 0 mod 25")
gg,hh = myHenselItere(f,g,h,s,t,m,2)
print(polynomeCentre(f-gg*hh,m^2)==0)
print("")
print("myHenselItere modulo 625")
print(myHenselItere(f,g,h,s,t,m,4))
gg,hh = myHenselItere(f,g,h,s,t,m,4)
print("A-t-in f-gg*hh = 0 mod 625")
print(polynomeCentre(f-gg*hh,m^2)==0)
print("")
print("Réponse à la question 5")
print(reponseQ5)

reset()
print("""\
# ****************************************************************************
# FACTORISATION AVEC LLL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p=13
k=4
m=p^k
j=3

PolZZ.<x> = PolynomialRing(ZZ)
f = x^4 - x^3 - 5*x^2 + 12*x - 6

u = x+7626

# Code pour l'EXERCICE

def polynomeCentre(f,m):
    Pol=f.parent()
    x=Pol.gen()
    coeff = f.list()
    f_tilde = Pol(0)
    for i in range(len(coeff)):
        c = coeff[i].mod(m)
        if c > m//2:
            c = c-m
        f_tilde =f_tilde + c*x^i
    retour = f_tilde
    return retour

def myHensel(f,g,h,s,t,m):
    Pol=f.parent()
    x=Pol.gen()
    e = polynomeCentre(f-g*h,m^2)
    q,r = (s*e).quo_rem(h)
    q,r = (polynomeCentre(q,m^2),polynomeCentre(r,m^2))
    g_etoile = polynomeCentre(g + t*e + q*g,m^2)
    h_etoile = h+r
    h_etoile = polynomeCentre(h_etoile,m^2)
    b = polynomeCentre(s*g_etoile +t*h_etoile-1,m^2)
    c,d = (s*b).quo_rem(h_etoile)
    c,d = polynomeCentre(c,m^2),polynomeCentre(d,m^2)
    s_etoile = polynomeCentre(s-d,m^2)
    t_etoile = polynomeCentre(t - t*b-c*g_etoile,m^2)
    retour = g_etoile,h_etoile,s_etoile,t_etoile
    return retour

def myHenselItere(f,g,h,m,l):
    Pol=f.parent()
    x=Pol.gen()
    d,ss,tt = xgcd(g,h)
    s=Pol(ss/mod(d,m));
    t=Pol(tt/mod(d,m))
    j = m
    ss= s
    tt = t
    ff = f
    gg=g
    hh = h
    while j < m^l:
        gg,hh,ss,tt = myHensel(ff,gg,hh,ss,tt,j)
        j = j^2
    retour = gg,hh
    return retour
    
# question 1a, on parcout simplement l'ensemble des éléments du corps
roots = []
for i in range(13):
    if f(i)%13 == 0:
        roots.append(i)
        
# on vérifie que roots a bien 4 éléments
alpha=roots[0]
beta=roots[1]
gamma=roots[2]
delta=roots[3]

# question 1b 
f1 = x-alpha
f2 = x-beta
f3 = x-gamma
f4 = x-delta

# on applique successivement le relèvement d'Hensel
gg1,ff4 = myHenselItere(f,f1*f2*f3,f4,p,k)
gg2,ff3 = myHenselItere(gg1,f1*f2,f3,p,k)
ff1,ff2 = myHenselItere(gg2,f1,f2,p,k)

alphahat=ff1.roots()[0][0]
betahat=ff2.roots()[0][0]
gammahat=ff3.roots()[0][0]
deltahat= ff4.roots()[0][0]

racine_mod = [alphahat,betahat,gammahat,deltahat]

# détermination de la base LLL
# le réseau est engendré par {u(x),u(x)x,m,mx,mx^2}

A = matrix(ZZ,3,5)
vec = [u,u*x,PolZZ(m),m*x,m*x^2]

for i in range(len(vec)):
    for j in range(3):
        A[j,i] = vec[i][j]

B = (A.T).LLL().T
base = []
v = vector([0,0,0])    #plus "petit" vecteur de la base LLL
for i in range(5):
    X = B[:,i]
    if X == 0:
        continue
    base.append(X)
    if v.norm()==0 or X.norm()<v.norm():
        v = X
    
# On peut alors envisage de prendre le plus petit vecteur de la base comme facteur
g = PolZZ(0)
v = v.list()
for i in range(len(v)):
    g = g + v[i]*x^i

#On peut remarquer que u est de degré 1, g est de degré 2 et f est de degré 4. On peut alors déduire simplement, le dernier facteur
h = f//(g*u)
factorisation =  Factorization([(u,1),(h,1),(g,1)])
# # Affichage des resultats

print("\n$1a/ Les racines sont", alpha, beta, gamma, delta,"modulo",p)
print("\n$1b/ Les racines sont", alphahat, betahat, gammahat, deltahat,"modulo",m)
print("")
print("Si on évalue f en ces valeurs, on a")
print("en alphahat ",f(alphahat))
print("en betahat",f(betahat))
print("en gammahat",f(gammahat))
print("en deltahat",f(deltahat))
print("")
print("u est-il un diviseur de f modulo 13^4 ?")
print(u.roots()[0][0] in racine_mod)
print("")
print("Une base LLL est ")
print(base)
print("")
print("On en déduit le facteur g suivant", g)
print("Est-il bien un facteur de f ?",f%g==0)
print("On obtient alors la factorisation suivante f=",factorisation)


