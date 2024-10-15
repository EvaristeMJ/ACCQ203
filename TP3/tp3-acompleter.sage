print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP3 : RESEAUX EUCLIDIENS                                                    #
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
# BASE D'UN RESEAU
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [ 2, -3,  4],
        [ 4,  4, 53]])

# Code pour l'EXERCICE
# L'ensemble est le noyau de la matrice (en retirant la dernière composante)
AA = matrix(ZZ,[
        [ 2, -3,  4,0],
        [ 1,  1, 2,-5]])
AA = AA.transpose()
L=[]
H,U = AA.hermite_form(transformation = True)
r = H.rank()
for i in range(r):
    L.append(U[i,:3])

# # Affichage des resultats
# Le réseau est bien un sous-groupe additif de (R^3,+) de la forme d'un réseau
print("\n$ Le réseau a pour base")
print(L)



reset()
print("""\
# ****************************************************************************
# APPLICATIONS NUMERIQUES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n1 = round(arctan(1),50)
n2 = round(arctan(1/5),50)
n3 = round(arctan(1/239),50)

r=-2.5468182768840820791359975088097915

Pol.<x>=PolynomialRing(ZZ)

# Code pour l'EXERCICE
M = 10^(6)
L1 = [1,0,0]
L2 = [0,1,0]
L3 = [0,0,1]
L4 = [round(M*n1),round(M*n2),round(M*n3)]
B = matrix(ZZ,[L1,L2,L3,L4]).T

alpha = B.LLL()[0]
r = -2.5468182768840820791359975088097915
r2 = r^2
r3 = r^3
M2 = 10^7
L1 = [1,0,0,0]
L2 = [0,1,0,0]
L3 = [0,0,1,0]
L4 = [0,0,0,1]
L5 = [round(M2),round(M2*r),round(M2*r2),round(M2*r3)]
B = matrix(ZZ,[L1,L2,L3,L4,L5]).T
beta = B.LLL()[0]
p = Pol(beta[3]*x^3+beta[2]*x^2+beta[1]*x+beta[0]) 
rr = p.roots()



# # Affichage des resultats

print("\n$ La relation de Machin est alpha1*n1+alpha2*n2+alpha3*n3=0 avec")
for i in range(3):
   print("alpha",i+1,"=",alpha[i])

print("\n$ Un polynome minimal plausible est")
print(p)
print("dont les racines sont")
print(p.roots(ring=RR,multiplicities=false))




reset()
print("""\
# ****************************************************************************
# ALGORITHME LLL
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

B = matrix(ZZ,[[  9, -25],
               [ -5,  14],
               [  4, -11]])

# Code pour l'EXERCICE

def myLLL(M):
    B = copy(M)
    n = B.ncols()

    repeter = true
    while repeter : 
        Betoile, mu = B.transpose().gram_schmidt(orthonormal=false)
        Betoile = Betoile.transpose()
        mu = mu.transpose()
        assert(B==Betoile*mu)

        for i in range(1,n):
            k = i-1
            while k>=0:
                B[:,i] = B[:,i]- round(mu[k,i])*B[:,k]
                mu[:,i] = mu[:,i]-round(mu[k,i])*mu[:,k]
                k=k-1
        repeter = False
        for i in range(n-1):
            if round(norm(Betoile[:,i])^2)>2*round(norm(Betoile[:,i+1])^2):
                B.swap_columns(i,i+1)
                repeter = True
                break
        assert(all(mu[i,j]<=1/2 for i in range(n) for j in range(i+1,n)))
        assert(all(mu[i,j]>=-1/2 for i in range(n) for j in range(i+1,n)))

    return B


# # Affichage des resultats

print("\n$ Une base LLL de B est")
print(myLLL(B))


reset()
print("""\
# ****************************************************************************
# RESEAUX CLASSIQUE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n = 8
e = vector([1/2]*8)

# Code pour l'EXERCICE

A = [1]*n
A = matrix(ZZ,A)
KernelA = A.transpose().kernel()
An = KernelA.matrix().transpose()
Ga = KernelA.gram_matrix()
an =Ga.determinant() #determinant de An

D = [1]*n+ [-2]
D = matrix(ZZ,D)
KernelD = D.transpose().kernel()
Dn = KernelD.basis_matrix()[:,:-1]
KernelD = Dn.row_module(base_ring = ZZ)
Dn = Dn.T
Gd = KernelD.gram_matrix()
dn = Gd.determinant()

E8 = matrix(QQ,Dn.nrows(),Dn.ncols()+1)
E8[:,:-1] = Dn
E8[:,-1] = matrix(QQ,e).T
ModE8 = E8.transpose().row_module(base_ring = ZZ)
e8 = ModE8.gram_matrix().determinant()
E8 = ModE8.basis_matrix().T


reponse6 = "Il suffit de regarder la base de E8. Elle montre qu'un vecteur de E8, s'écrit n*e + x, où x est un vecteur d'entiers. Si n est pair, n*e est un vecteur d'entiers et donc n*e+x aussi. Si n est impair, n*e est un vecteur de demi entiers, et donc n*e + x l'est aussi."
# Vecteurs minimaux de An
# Pour avoir norm(x) = 2, on doit avoir au plus deux coefficients non nuls et égaux à 1. Puisque la somme des coefficients doit être nulles, on en déduit aussi que le vecteur est de la forme (0,...,1,0,...,0,-1,0...)
reponse7_an = []
for i in range(n):
    for j in range(n):
        if i != j:
            x = vector(ZZ,n)
            x[i] = 1
            x[j] = -1
            
            reponse7_an.append(x)
            
# Vecteurs minimaux de Dn
# Pour avoir norm(x) = 2, on doit avoir au plus deux coefficients non nuls. Puisque la somme des coefficients doit êtres nulles, on en déduit aussi que le vecteur est de la forme (0,...,+-1,0,...,0,+-1,0...) (donc les vecteurs minimaux de An sont aussi des vecteurs minimaux de Dn)
reponse7_dn = copy(reponse7_an)
for i in range(n):
    for j in range(i,n):
        if j!=i:
            x = vector(ZZ,n)
            x[i] = 1
            x[j] = 1
            reponse7_dn.append(x)
            x = vector(ZZ,n)
            x[i] = -1
            x[j] = -1
            reponse7_dn.append(x)
            
#Vecteurs minimaux de E8
# Si les vecteurs sont à coordonnées entières, idem que pour Dn, si jamais les coordonnées ne sont pas entières les coordonnées ne peuvent être égales qu'à +- 1/2
# Dans ce cas, on a pour un vecteur coordonnées +- 1/2, une norme de 2
# Désormais, il suffit que la somme des coordonnées soit paire. Pour cela, il faut que les moins aillent par pair.

reponse7_e8 = copy(reponse7_dn)

reponse7_e8.append(e)
reponse7_e8.append(-e)


for i in range(n):
    for j in range(i,n):
        if j!=i:
            x = vector(QQ,n)
            x = x + e
            x[i] = -1/2
            x[j] = -1/2
            reponse7_e8.append(x)
            reponse7_e8.append(-x)

for i in range(n):
    for j in range(i,n):
        for k in range(j,n):
            for l in range(k,n):
                if l!=k and k!= j and j!= i:
                    x = vector(QQ,n)
                    x = x + e
                    x[i] = -1/2
                    x[j] = -1/2
                    x[k] = -1/2
                    x[l] = -1/2
                    reponse7_e8.append(x)

# Les vecteurs minimaux sont ceux de norme égale à 2. Il n'existe pas d'élément de norme strictement inférieure à 2 dans An, Dn ainsi que E8. Puisque les normes sont entières, la seule valeur possible serait 1. Hors dans les trois cas, il n'existe pas d'éléments de norme unitaire.
# On représente le nombre de vecteurs minimaux des trois réseaux sous la forme (An,Dn,E8)
# On a pour An : n(n-1) Dn : 2n(n-1) et E8 : 240
reponse8 = (len(reponse7_an),len(reponse7_dn),len(reponse7_e8))

# # Affichage des resultats

print("\n$ nombre d'éléments minimaux pour An, Dn et E8")
print(reponse8)

print("\n$ Une base de An est")
print(An, "de déterminant",an)

print("\n$ Une base de Dn est")
print(Dn, "de déterminant",dn)

print("\n$ Une base de E8 est")
print(E8, "de déterminant",e8)


reset()
print("""\
# ****************************************************************************
# DENSITES OPTIMALES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE
A = [1]*2
A = matrix(ZZ,A)
A2 = A.transpose().kernel()
discA2 = sqrt(A2.gram_matrix().determinant())

B = [1]*3
B = matrix(ZZ,B)
A3 = B.transpose().kernel()
discA3 = sqrt(A3.gram_matrix().determinant())

D = [1]*4+ [-2]
D = matrix(ZZ,D)
KernelD = D.transpose().kernel()
D = KernelD.basis_matrix()[:,:-1]
D4 = D.row_module(base_ring = ZZ)
discD4 = sqrt(D4.gram_matrix().determinant())

D = [1]*5+ [-2]
D = matrix(ZZ,D)
KernelD = D.transpose().kernel()
D = KernelD.basis_matrix()[:,:-1]
D5 = D.row_module(base_ring = ZZ)
discD5 = sqrt(D5.gram_matrix().determinant())

e = vector([1/2]*8)

D = [1]*8+ [-2]
D = matrix(ZZ,D)
KernelD = D.transpose().kernel()
D = KernelD.basis_matrix()[:,:-1]
D8 = D.row_module(base_ring = ZZ)

E = matrix(QQ,D8.matrix().ncols()+1,D8.matrix().nrows())
E[:-1,:] = D8.matrix()
E[-1,:] = matrix(QQ,e)
E8 = E.row_module(base_ring = ZZ)
discE8 = sqrt(E8.gram_matrix().determinant())

def densite(lambda_min,disc,n):
    return sqrt((pi*lambda_min)^n)/(disc*factorial(n/2)*2^n)

lambda_min = 2

a2 = densite(lambda_min,discA2,2)
a3 = densite(lambda_min,discA3,3)
d4 = densite(lambda_min,discD4,4)
d5 = densite(lambda_min,discD5,5)
e8 = densite(lambda_min,discE8,8)

# # Affichage des resultats

print("\n$ La densité de A2 est",a2)
print("\n$ La densité de A3 est",a3)
print("\n$ La densité de D4 est",d4)
print("\n$ La densité de D5 est",d5)
print("\n$ La densité de E8 est",e8)




