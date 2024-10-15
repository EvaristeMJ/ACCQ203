print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP11 : CODES CORRECTEURS D'ERREURS ALGEBRIQUES                              #
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
# CODE DE HAMMING
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

C = codes.HammingCode(GF(2),3)

# Code pour l'EXERCICE

B = C.generator_matrix() # a completer
d = C.minimum_distance() # a completer

# # Affichage des resultats

print( "Le code",C,"a pour matrice")
print( B)
print( "et pour distance minimale d=",d)


reset()
print("""\
# ****************************************************************************
# DECODAGE DE BERLEKAMP-MASSEY
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

from sage.matrix.berlekamp_massey import berlekamp_massey

q=13
Fq = FiniteField(q)
k=5
alpha=Fq(6)
PolFq.<x> = PolynomialRing(Fq)
r = vector(Fq,[5,2,9,9,1,11,7,5,7,4,3,1])

# Code pour l'EXERCICE

G = matrix()
H = matrix()
s = vector(Fq,0)
sigma = PolFq(0)
M = []
m = vector(Fq,0)

# # Affichage des resultats

print( "La matrice génératrice du code est")
print( G)
print( "La matrice de controle du code est")
print( H)
print( "Le syndrome est")
print( s)
print( "Le polynome localisateur d'erreurs est")
print( sigma)
print( "La position des erreurs est")
print( M)

reset()

print("""\
# ****************************************************************************
# DECODAGE DE BERLEKAMP-MASSEY
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

from sage.matrix.berlekamp_massey import berlekamp_massey

q=13
n=q-1
Fq = FiniteField(q)
k=5
alpha=Fq(6)
PolFq.<x> = PolynomialRing(Fq)
r = vector(Fq,[5,2,9,9,1,11,7,5,7,4,3,1])

# Code pour l'EXERCICE

def G_matrix(X,k,n):
    G = matrix(k,n)
    for i in range(k):
        for j in range(n):
            G[i,j] = X[j]^i
    return G

def H_matrix(X,k,n):
    H = matrix(n-k+1,n)
    for i in range(n-k+1):
        for j in range(n):
            H[i,j] = X[j]^i
    return H

X = [alpha^i for i in range(0,n)]  
            
G = G_matrix(X,k,n)
H = H_matrix(X,k,n)
s = H*r
sigma = berlekamp_massey(list(s)).reverse()
M = []
for i in range(n):
    if sigma(alpha^(-i))==Fq(0):
        M.append(i)
        
GG = G.matrix_from_columns([i for i in range(n) if i not in M ])
rr = vector([r[i] for i in range(len(r)) if i not in M])

print(GG)
print(rr)
m = GG.solve_left(rr)

# # Affichage des resultats

print( "La matrice génératrice du code est")
print( G)
print( "La matrice de controle du code est")
print( H)
print( "Le syndrome est")
print( s)
print( "Le polynome localisateur d'erreurs est")
print( sigma)
print( "La position des erreurs est")
print( M)
print( "Le message envoye le plus probable est")
print( m)



print( "Le message envoye le plus probable est")
print( m)



reset()
print("""\
# ****************************************************************************
# DECODAGE EN LISTE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=23
Fq = FiniteField(q)
k = 3
alpha=Fq(14)
MPol.<x,y>=PolynomialRing(Fq,2)
n = q-1
r=vector(Fq, [12,18,15,22,17,5,14,21,17,4,13,8,4,10, 15,11,22,12,13,9,14,12])

# Code pour l'EXERCICE

def to_char(m):
    return [chr(Integer(m[i])+ord('a')) for i in range(len(m))]

def hamming(x,y):
    return len([i for i in range(len(x)) if x[i]!= y[i]])

b1= floor(((q - k) - 1)/2)
b2=floor(n-sqrt(2*(k-1)*n))

exposants = []
for i in range(b2):
    for j in range(b2):
        if i+2*j < b2:
            exposants.append((x^i) * (y^j))
            
M = matrix(len(exposants),len(r))
for i in range(len(r)):
    for j in range(len(exposants)):
        f = exposants[j]
        M[j,i] = f(alpha^i,r[i])

K = kernel(M)
pol_basis = []

for p in K.basis() :
    pol_basis.append(sum([p[i]*exposants[i] for i in range(len(p))]))

facteur = []    


for p in pol_basis :
    for pp,_ in p.factor():
        if pp.polynomial(y).degree() == 1 :
            fact = (-pp+y).polynomial(x)
            card = len([i for i in range(len(r)) if fact(alpha^i) == r[i]]) #indices où on a égalité de fact(alpha^i) et r_i
            
            if fact.degree() < k and card  >= q-1-b2 :
                facteur.append(fact)

facteur = set(facteur)

L = [vector(Fq,[p(alpha^i) for i in range(q-1)]) for p in facteur]
G = matrix([[z^j for z in [alpha^i for i in range(22)]] for j in range(k)]) 
M = [to_char(G.solve_left(c)) for c in L]
c = vector(Fq,0)
m = ""

if hamming(L[0],r) >= hamming(L[1],r):
    c = L[0]
    m = M[0]
else:
    c = L[1]
    m = M[1]



# # Affichage des resultats

print( "Le decodage unique fonctionne jusqu'à distance",b1)
print( "Le decodage de Sudan fonctionne jusqu'à distance",b2)
print( "La liste des mots de code les plus proches de r est",L)
print( "Elle correspond aux messages")
print( M)
print( "Le mot de code le plus probable est",c)
print( "soit le message",m)

reset()

print("""\
# ****************************************************************************
# BORNES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=49

# Code pour l'EXERCICE

def H(d):
    if d == 0:
        return 0
    return d*log(q-1,q)-d*log(d,q) - (1-d)*log(1-d,q)
    
def gv(d):
    if d<= 1-1/q:
        return 1-H(d)
    return 0

def plotkin(d):
    theta = 1-1/q
    if d <= theta:
        return 1-d/theta
    return 0

def hamming(d):
    return 1-H(d/2)

def mrrw(d):
    return H((q-1-(q-2)*d-2*sqrt((q-1)*d*(1-d)))/q)

# # Affichage des resultats

#plot(D,gv,xmin=0,xmax=1,ymin=0,ymax=1)
plot([gv,plotkin,hamming,mrrw],xmin=0,xmax=1,ymin=0,ymax=1, legend_label=["GV","Plotkin","Hamming","MRRW"],title="Differentes bornes sur les codes")



print("""\
# ****************************************************************************
# BORNES BATTANT GILBERT-VARSHAMOV
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE

def codesalg(d):
    return 1-d-1/(sqrt(q-1))

# # Affichage des resultats

plot([gv,plotkin,hamming,mrrw,codesalg],xmin=0,xmax=1,ymin=0,ymax=1, legend_label=["GV","Plotkin","Hamming","MRRW","Codes Algebriques"],title="Differentes bornes sur les codes")
