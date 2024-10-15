print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP14 : LOG DISCRET ET COUPLAGES                                             #
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
# PAS DE BEBE, PAS DE GEANT
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p1 = 1823
Fp1 = FiniteField(p1)
b1 = Fp1(3)
x1 = Fp1(693)

p2 = 239
Fp2 = FiniteField(p2)
b2 = Fp2(2)
x2 = Fp2(15)


# Code pour l'EXERCICE

def Shanks(x,b):
    Fp = x.parent()
    p = Fp.cardinality()
    s = floor(sqrt(p-1))+1
    T = [-1]*(p-1)
    for j in range(s):
        beta = x*b^(-j)
        T[beta] = j
    i = 0
    gamma = 1
    bb = b^s
    while True:
        i = i +1
        gamma = gamma*bb
        if T[gamma]>= 0:
            break
    return i*s+T[gamma]


# # Affichage des resultats

print("Question 2 :", Shanks(x1,b1))
print("Question 3 :", Shanks(x2,b2))




reset()
from random import randint
print("""\
# ****************************************************************************
# RHO DE POLLARD
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p= 281
Fp = FiniteField(p)
x1 = Fp(263)
b1 = Fp(239)
x2 = Fp(165)
b2 = Fp(127)
x3 = Fp(210)
b3 = Fp(199)


# Code pour l'EXERCICE

def rho(g,b):
    Fp = g.parent()
    p = Fp.cardinality()
    n = b.multiplicative_order()
    Zn = Integers(n)
    def phi(alpha,beta,w):
        h = hash(w)
        h = h%3
        if h==0:
            return (alpha,beta+1,g*w)
        elif h== 1:
            return (2*alpha,2*beta,w^2)
        else:
            return (alpha+1,beta,b*w)
    
    while  True:
        alpha,beta = ZZ.random_element(1,n),ZZ.random_element(1,n)
        ax,bx,x = phi(alpha,beta,(b^alpha)*(g^beta))
        ay, by, y = phi(ax,bx,x)

        while x != y:
            ax,bx,x = phi(ax,bx,x)

            ay, by, y = phi(ay,by,y)
            ay, by, y = phi(ay,by,y)
        try:
            return Zn((ax-ay).mod(n)*inverse_mod(by-bx,n))
        except:
            continue


# # Affichage des resultats

print("Le log de x=",x1,"en base",b1,"vaut",rho(x1,b1),".")
print("Le log de x=",x2,"en base",b2,"vaut",rho(x2,b2),".")
print("Le log de x=",x3,"en base",b3,"vaut",rho(x3,b3),".")



reset()
print("""\
# ****************************************************************************
# COUPLAGE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[11,0])
print("groupe de E=", E.abelian_group()) # pour verifier
S = E(24,34)
T = E(5,27)
r = 10
print("Verification de la r-torsion : r*S =", r*S, "et r*T =", r*T)


# Code pour l'EXERCICE

def myLine(P1,P2,S):
    E=P1.curve()
    K = E.base_field()
    alpha = E.a4() # alpha tel que y^2 = x^3 + alpha x +beta (equation de weierstrass)
    x1=P1[0]; y1=P1[1]; z1=P1[2]
    x2=P2[0]; y2=P2[1]; z2=P2[2]
    xS=S[0]; yS=S[1]; zS=S[2]
    a = 1
    
    # distinction de cas pour le calcul de g(P1,P2)(S) qui représente les équations de droite
    if z1 == 0 and z2 == 0:
        return K(1)
    
    elif z1 == 0:
        return xS-x2*zS
    elif z2 == 0 or (x1==x2 and y1 == -y2):
        return xS-x1*zS

    elif x1==x2:
        a = (3*x1^2+alpha)/(2*y1)
    
    else:
        a = (y1-y2)/(x1-x2)
    return yS-y1*zS - (xS-x1*zS)*a


def myH(P1,P2,S):
    return myLine(P1,P2,S)/(myLine(P1+P2,-P1-P2,S))

def myMiller(r,S,P):
    R = S
    br = r.bits()
    l = len(br)
    f = 1
    for i in range(l-2,-1,-1):
        f = myH(R,R,P)*(f^2)
        R = 2*R
        if br[i] == 1:
            f *= myH(R,S,P)
            R += S
    return f

def myTatePairing(S,T,r):
    E = S.curve()
    try:
        return myMiller(r,S,T)
    except ZeroDivisionError:
        Q = E.random_point()
        return myTatePairing(S,T+Q,r)


# # Affichage des resultats
print("Calcul du couplage", myTatePairing(S,T,r))
print("D'après Sage, le couplage est égale à :", S.tate_pairing(T,r,1))



reset()
print("""\
# ****************************************************************************
# ATTAQUE M.O.V.
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2199023255579
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,0])
P = E(1435967701832 , 123951463462)
Q = E(1129476910351 , 1383670460733)

# Code pour l'EXERCICE

j = E.j_invariant() # j-invariant a faire calculer par une fonction de SageMath
rep2 = f"Le j-invariant est égal à {j}et p = {p} est congru à 3 modulo 4 donc la courbe est supersingulière."
t = 1 # Ecrire le code pour calculer cette valeur
r = P.order()

# Calcul de t
while (p^t-1)%r != 0:
    t = t+1

q = p^t
Fq.<alpha> = FiniteField(q)
EE = EllipticCurve(Fq,[1,0])
PP = EE(1435967701832 , 123951463462)
QQ = EE(1129476910351 , 1383670460733)
SS = EE.random_element() # point a calculer vous-meme
rr = PP.order()
while SS.order() != rr or PP.weil_pairing(SS,r) == 1:
    SS = EE.random_element()
zeta1 = Fq(PP.weil_pairing(SS,r))
zeta2 = Fq(QQ.weil_pairing(SS,r))
lambd = log(zeta2,zeta1)


# # Affichage des resultats

print("p premier ?",p.is_prime())
print("j-invariant de E :",j)
print("p mod 4 =", mod(p,4))
print(rep2)
print("Cardinal de E(Fp) :",E.cardinality(),"=",E.cardinality().factor())
print("Ordre de P :",P.order())
print("Cardinal de E(Fq) :",EE.cardinality(),"=",EE.cardinality().factor())
print("Point S :",SS)
print("On calcule zeta1 =",zeta1,", zeta2 =",zeta2,", lambda =",lambd,".")




reset()
print("""\
# ****************************************************************************
# CALCUL D'INDICE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 439
Fp = FiniteField(p)
g = Fp(237)
b = Fp(136)
y = 11

# Code pour l'EXERCICE

def factorisation_friable(x,P):
    """
    Si x est y-friable, renvoie (True, exposants de la factorisation dans la base de friabilité)
    Sinon renvoie (False, None)
    """
    if x==0:
        return (False,None)
    factorisation = ZZ(x).factor()
    ex = vector([0]*len(P))
    for (f,e) in factorisation:
        if f>P[-1]:
            return (False,None)
        ex[P.index(f)] = e
    return (True,ex)

def LogIndice(g,b,y):
    Fp = g.parent()
    p = Fp.cardinality()
    P = [i for i in range(2,y+1) if ZZ(i).is_prime()]
    k = len(P)
    Zn = Integers(p-1)
    M = (ZZ^k)/((p-1)*ZZ^k)

    A = []
    V = []

    while True:
        A = []
        V = []
        while len(V) <= 4*k:
            alpha = ZZ.random_element(0,p-1)
            gamma = ZZ((b^alpha).mod(p))
            
            res = factorisation_friable(gamma,P)
            if res[0]:
                A.append(alpha)
                V.append(res[1])
        if M.submodule(V) == M:
            break

    #Résolution du système pour calculer les logarithmes
    B = matrix(Zn,V)
    a = vector(Zn,A)
    L = B.solve_right(a)
    
    while True:
        beta = ZZ.random_element(0,p-1)
        res = factorisation_friable(g*b^beta,P)
        if res[0]: #On s'arrête quand g*b^beta est y-friable
            break
    expo = res[1]
    l = -beta
    for i in range(len(P)):
        if expo[i] == 0:
            continue
        l += expo[i]*L[i]
    return l

# # Affichage des resultats

print("Le log de g=",g,"en base",b,"vaut",LogIndice(g,b,y),".")
print("Et le vrai résultat est : ", log(g,b))


