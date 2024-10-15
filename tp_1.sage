
print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP1 : INTRODUCTION A SAGEMATH                                               #
# *************************************************************************** #
# *************************************************************************** #
""")

# CONSIGNES
#
# Les seules lignes a modifier sont annoncee par "Code pour l'exercice n"
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
# EXERCICE 1
# ****************************************************************************
""")

# Code pour l'exercice 1

p = random_prime(10000)
q = random_prime(10000)
n = p*q
phi = (p-1)*(q-1)

e1 = ZZ.random_element(phi)


while (true):
    try:
        e2 = ZZ.random_element(phi)
        e2 = mod(e2,phi)
        d=1/e2
        break
    except ZeroDivisionError:
        pass


def chiffrement(a):
    return mod(a,n)^e2
dechiffrement = lambda b : mod(b,n)^d

for a in range(4):
    dechiffrement(chiffrement(a))

reponse = "Votre reponse ici." # A MODIFIER

# Affichage des resultats pour l'exercice 1

print("On a choisi p=",p,"q=",q)
print("q est-il premier ?", q.is_prime())
if p==q:
    print('Mauvais choix')
else:
    print('Choix correct')
print("Type de e1 : ",type(e1))
print("Parent de e1 : ",parent(e1))
print("Inverse de e1 :",1/e1)
print("Type de e2 : ",type(e2))
print("Parent de e2 : ",parent(e2))
print("Inverse de e2 :",1/e2)

print("\n$ Question")
print(reponse)
# algorithm RSA


reset()
print("""
# ****************************************************************************
# EXERCICE 2
# ****************************************************************************
""")

# Code pour l'exercice 2

a = 119
b = 435
d,u,v = xgcd(a,b)

alpha = 2
beta = 3
m = 45
n = 14

x = crt(alpha,beta,m,n) # A MODIFIER

# Affichage des resultats pour l'exercice 2
print("\n$ Question 1")
print("Une relation de Bezout entre a=",a," et b=",b," est au + bv = d avec")
print("d=",d)
print("u=",u)
print("v=",v,".")

print("\n$ Question 2")
print("Une solution au systeme")
print("x=",alpha,"mod",m)
print("x=",beta,"mod",n)
print("est x=",x,".")

reset()
print("""
# ****************************************************************************
# EXERCICE 3
# ****************************************************************************
""")

# Code pour l'exercice 3

p = 157
k = 4
q = p^k
Fp = FiniteField(p)
Fq.<alpha> = FiniteField(q)
z = 130*alpha^3 + 97*alpha^2 + 99*alpha + 18

vectors = []

for i in range(k):
	vectors.append(vector(alpha^(i*p))

# Construction de Q A COMPLETER
Q = matrix(Fp,vectors).T


t = vector([1,alpha,alpha^2,alpha^3])

sigmaz = z^p
sigmaz = t*Q*vector(z)

test1 = True
for _ in range(10):
	x = Fq.random_element()
	sigmax1 = x^p
	sigmax2 = t*Q*vector(x)
	if sigmax1 != sigmax2:
		test1 = False
		break
print(test1)

def tau1(x):
    """
    Premiere facon de calculer la racine p-ieme.
    """
    x=Fq(x)
    return x.pth_root(p)

def tau2(x):
    """
    Seconde facon de calculer la racine p-ieme.
    """
    x=Fq(x)
    return t*Q.inverse()*vector(x)

tau1z = tau1(z)
tau2z = tau2(z)

test2 = True
for _ in range(100):
	r = Fq.random_element()
	if tau1(r) != tau2(r):
		test2 = False
		break


# Affichage des resultats pour l'exercice 3

print("$ Question 1")
print("La matrice du Frobenius vaut Q=")
print(Q)
print("La valeur de sigma(z) est",sigmaz)
print("Verification de Q sur une centaine d'exemple :",test1)
print("$ Question 2")
print("La valeur de tau(z) est")
print(tau1z)
print("ou")
print(tau2z)
print("Verification de tau1 et tau2 sur une centaine d'exemple :",test2)


reset()
print("""
# ****************************************************************************
# EXERCICE 4
# ****************************************************************************
""")

# Code pour l'exercice 4
# Question 1

p = random_prime(8)
p = 2
n = ZZ(randint(2,4)) # randint renvoie un entier de python
n = 4
Fp = FiniteField(p)
PolFp.<x> = PolynomialRing(Fp)
t = p^n - 1

L_irred = []

for P in PolFp.polynomials( of_degree = n):
	if P.is_irreducible():
		L_irred.append(P)

# Question 2

L_prim = []
for P in L_irred:
	c = P.coefficients()
	if c[0] == 1 and P.is_primitive():
		L_prim.append(P)


# Question 3

l_prim = len(L_prim)
test1 = False
if L_prim == euler_phi(t)/n:
	test1 = True

Phi = cyclotomic_polynomial(t)
P = L_prim[0]*L_prim[1]
test2 = False
if P == Phi:
	test2 = True

# Question 4


#p1 = L_irred[0]
#p2 = L_irred[1]
reponse_4a="A ECRIRE"
reponse_4b="A ECRIRE"
#Fq1.<alpha1> = FiniteField(p^n,modulus=p1)
#Fq2.<alpha2> = FiniteField(p^n,modulus=p2)
Psi = 0
#Psi = Fq1.hom([0])
psi  = 0
#psi  = Psi(alpha1^2)


# Affichage des resultats pour l'exercice 4

print("\n$ Question 1")
print("Les polynomes irreductibles unitaires de degre",n,"de", PolFp,"sont les suivants :")
for p in L_irred:
    print(p)

print("\n$ Question 2")
print("Les polynomes primitifs unitaires de degre",n,"de", PolFp,"sont les suivants :")
for p in L_prim:
    print(p)

print("\n$ Question 3a")
print("Il y a ",l_prim," polynomes primitifs.")
print("Ce compte est-il egal a phi(t)/n ?",test1)
print("$ Question 3b")
print("Verifie-t-on la formule du produit ?",test2)

print("\n$ Question 4a")
print(reponse_4a)
print("$ Question 4b")
print(reponse_4b)
print("$ Question 4c")
print("Un morphisme possible est")
print(Psi)
print("qui envoie alpha1^2 sur")
print(psi)

reset()
print("""
# ****************************************************************************
# EXERCICE 5
# ****************************************************************************
""")

# Donnees de l'exercice 5

A = matrix(QQ,[
        [-2, 1,  1],
        [ 8, 1, -5],
        [ 4, 3, -3]])
C = matrix(QQ,[
        [  1,  2, -1],
        [  2, -1, -1],
        [ -5,  0,  3]])

# Code pour l'exercice 5

W = C.right_kernel()
test=True
# TEST A ECRIRE

X0 = 0 # A MODIFIER



reponse_3 = "A ECRIRE"
# Affichage des resultats pour l'exercice 5

print("\n$ Question 1")
print("ker C est inclus dans ker A :", test)
print("\n$ Question 2")
print("Solution particuliere :")
print(X0)
print("\n$ Question 3")
print("Solution generale :", reponse_3)

reset()
print("""
# ****************************************************************************
# EXERCICE 7
# ****************************************************************************
""")

# Donnees de l'exercice 7

A = matrix(ZZ,3,3,[
        [ 3, -7,  2],
        [ 4,  6, -1],
        [-4,  3,  2]])

# Code pour l'exercice 7

print("forme LU") # A MODIFIER
print("determinant de A") # A MODIFIER
reponse = "A ECRIRE"

# Affichage des resultats pour l'exercice 7

print(reponse)
