print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP2 : ALGEBRE LINEAIRE SUR UN ANNEAU PRINCIPAL                              #
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
# MISE SOUS FORME NORMALE D'HERMITE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [-2,  3,  3,  1],
        [ 2, -1,  1, -3],
        [-4,  0, -1, -4]])

A1 = random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=3)

U = identity_matrix(4)

# Code pour l'EXERCICE

def cherche_pivot_non_nul(i,j,A,U):
    """
    Echange la colonne j avec une colonne à gauche pour que A[i,j] soit non nul.
    """
    for k in range(0,j):
        if A[i,k]!=0:
            A.swap_columns(k,j)
            U.swap_columns(k,j)
            break

def normalise_pivot(i,j,A,U):
    """
    Multiplie la colonne j si besoin pour que A[i,j] soit positif.
    """
    if A[i,j]<0:
        A[:,j] = -A[:,j]
        U[:,j] = -U[:,j]
    assert(A[i,j]>0)

def annule_a_gauche(i,j,A,U):
    """
    Annule les coefficients à gauche de A[i,j]
    """
    for k in range(j):
        (d,u,v)=xgcd(A[i,j],A[i,k])
        s = -A[i,k]/d
        t = A[i,j]/d
        Cj = A[:,j]
        Uj = U[:,j]
        Ck = A[:,k]
        Uk = U[:,k]
        A[:,k] = s*Cj+t*Ck
        U[:,k] = s*Uj+t*Uk
        A[:,j] = u*Cj + v*Ck
        U[:,j] = u*Uj + v*Uk

def reduit_a_droite(i,j,A,U):
    """
    Réduit les coefficients à gauche de A[i,j] modulo A[i,j]
    """
    n = A.nrows()
    for k in range(j+1,n+1):
        q = A[i,k]//A[i,j]
        A[:,k] = A[:,k]-q*A[:,j]
        U[:,k] = U[:,k] - q*U[:,j]


def MyHNF(A):
    """
    Forme normale d'Hermite selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    U = identity_matrix(n)
    # ECRIVEZ VOTRE CODE ICI, VOUS POUVEZ REPRENDRE LES FONCTIONS PRECEDENTES
    # COMME SOUS-FONCTION
    l = max(1,m-n+1)-1
    i = m-1
    j = n-1
    while i >= l:
        cherche_pivot_non_nul(i,j,H,U)
        if H[i,j] == 0:
            i = i - 1
        else:
            normalise_pivot(i,j,H,U)
            while H[i,:j] != 0:
                normalise_pivot(i,j,H,U)
                annule_a_gauche(i,j,H,U)
            reduit_a_droite(i,j,H,U)
            i = i-1
            j = j-1
                
    assert(H-A*U==0)
    return H,U

def SageHNF(A):
    """
    Forme normale d'Hermite de SAGE avec la convention francaise :
    Les vecteurs sont ecrits en colonne.
    """
    m = A.nrows()
    n = A.ncols()
    Mm = identity_matrix(ZZ,m)[::-1,:]
    Mn = identity_matrix(ZZ,n)[::-1,:]
    AA = (Mm*A).transpose()
    HH, UU = AA.hermite_form(transformation=True)
    H = (HH*Mm).transpose()*Mn
    U = UU.transpose()*Mn
    assert(H-A*U==0)
    return H,U

H,  U  = MyHNF(A)
HH, UU = SageHNF(A)

test = True # Test a ecrire

for _ in range(100):
    A1 = random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=3)
    H1,U1 = MyHNF(A1)
    HH1,UU1 = SageHNF(A1)
    if H1-HH1 != 0:
        test = false
        break

# # Affichage des resultats

print("\n$ Question 4")
print("La matrice A = ")
print(A)
print("a pour forme normale d'Hermite H=")
print(H)
print("et matrice de transformation U=")
print(U)
print("\n$ Question 5")
print("D'apres SageMath, la matrice A a pour forme normale d'Hermite H=")
print(HH)
print("et matrice de transformation U=")
print(UU)
print("\n$ Question 6")
print("Les deux fonctions coincident-elles sur une centaine d'exemples ?")
print(test)
reset()
print("""\
# ****************************************************************************
# MISE SOUS FORME NORMALE DE SMITH
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,2,3,[
        [4, 7, 2],
        [2, 4, 6]])

X2 = matrix(ZZ,3,3,[
        [-397, 423, 352],
        [   2,  -3,   1],
        [-146, 156, 128],
])

PolQ.<xQ> = PolynomialRing(QQ)
AQ = matrix(PolQ,3,[
            [xQ + 1,  2,     -6],
            [     1, xQ,     -3],
            [     1,  1, xQ - 4]])
Pol2.<x2> = PolynomialRing(FiniteField(2))
AF2 = matrix(Pol2,3,[
            [x2 + 1,  2,     -6],
            [     1, x2,     -3],
            [     1,  1, x2 - 4]])
Pol3.<x3> = PolynomialRing(FiniteField(3))
AF3 = matrix(Pol3,3,[
            [x3 + 1,  2,     -6],
            [     1, x3,     -3],
            [     1,  1, x3 - 4]])
Pol5.<x5> = PolynomialRing(FiniteField(5))
AF5 = matrix(Pol5,3,[
            [x5 + 1,  2,     -6],
            [     1, x5,     -3],
            [     1,  1, x5 - 4]])

# Code pour l'EXERCICE

def chercher_pivot(k,A,L,R):
    m = A.nrows()
    n = A.ncols()
    if A[k,k] != 0:
        return
    else:
        for i0 in range(k,m):
            for j0 in range(k,n):
                if A[i0,j0] != 0:
                    A.swap_lines(k,i0)
                    L.swap_lines(k,i0)
                    A.swap_columns(k,j0)
                    R.swap_columns(k,j0)

def annuler_en_bas(k,i,A,L,R):
    (d,s,t) = xgcd(A[k,k],A[i,k])
    if A[i,k]%A[k,k]==0:
        d = A[k,k]
        s = 1
        t = 0
    u = -A[i,k]/d
    v = A[k,k]/d
    Lk = A[k,:]
    LLk = L[k,:]
    Li = A[i,:]
    LLi = L[i,:]
    A[k,:] = s*Lk + t*Li
    L[k,:] = s*LLk + t*LLi
    A[i,:] = u*Lk + v*Li
    L[i,:] = u*LLk + v*LLi

def annuler_a_droite(k,j,A,L,R):
    (d,s,t) = xgcd(A[k,k],A[k,j])
    if A[k,j]%A[k,k]==0:
        d = A[k,k]
        s = 1
        t = 0
    u = -A[k,j]/d
    v = A[k,k]/d
    Ck = A[:,k]
    CRk = R[:,k]
    Cj = A[:,j]
    CRj = R[:,j]
    A[:,k] = s*Ck + t*Cj
    R[:,k] = s*CRk + t*CRj
    A[:,j] = u*Ck + v*Cj
    R[:,j] = u*CRk + v*CRj

def divise_le_reste(k,A,L,R):
    m = A.nrows()
    n = A.ncols()
    for j in range(k+1,n):
        for i in range(k+1,m):
            if A[i,j]%A[k,k]!=0:
                A[:,k] = A[:,k] + A[:,j]
                R[:,k]= R[:,k] + R[:,j]
                return

def coefficients_non_nuls(k,A):
    m = A.nrows()
    n = A.ncols()
    for j in range(k+1,n):
        if A[k,j] != 0:
            return True
    for i in range(k+1,m):
        if A[i,k] != 0:
            return True
    return False

def MySNF(A):
    """
    Forme normale de Smith selon votre code
    """
    t = A.base_ring() #Permet de définir une matrice identité dans l'anneau considéré
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    L = identity_matrix(t,m)
    U = identity_matrix(t,n)

    for k in range(min(m,n)):
        while coefficients_non_nuls(k,H):
            chercher_pivot(k,H,L,U)
            for i in range(k+1,m):
                annuler_en_bas(k,i,H,L,U)
            for j in range(k+1,n):
                annuler_a_droite(k,j,H,L,U)
            divise_le_reste(k,H,L,U)
    assert(H-L*A*U==0)
    return H,L,U

H1, L1, U1 = MySNF(X1)
H2, L2, U2 = MySNF(X2)

HQ, _, _ = MySNF(AQ)
HF2, _, _ = MySNF(AF2)
HF3, _, _ = MySNF(AF3)
HF5, _, _ = MySNF(AF5)

smith_forms = [H1,H2,HQ,HF2,HF3,HF5]
matrices = [X1,X2,AQ,AF2,AF3,AF5]

def sont_associes(a,b):
    """
    Puisqu'on travaille dans des anneaux principaux, il suffit de montrer que a divise b et b divise a
    Renvoie vrai si a et b sont associés
    """
    if b==0 and a==0:
        return True
    
    if b!= 0 and a%b==0:
        if a!=0 and b%a == 0:
            return True
    return False
    
def meme_SNF(H1,H2):
    """
    Renvoie vrai si les coefficients de la SNF sont tous associés
    Et donc que les coefficients de la SNF sont identiques à un facteur inversible près
    Renvoie faux si les deux SNF ne sont pas équivalentes
    """
    if H1.nrows() != H2.nrows():
        return False
    if H1.ncols() !=  H1.ncols():
        return False
    n = H1.ncols()
    m = H1.nrows()
    for k in range(min(n,m)):
        if not sont_associes(H1[k,k],H2[k,k]):
            return False
    return True

test = True
for i in range(len(matrices)):
    A = matrices[i]
    H = smith_forms[i]
    smith_sage = A.smith_form()[0]
    if not meme_SNF(H,smith_sage):
        test=False
        break

# # Affichage des resultats

print("\n$ Question 4")
print("La matrice X1 = ")
print(X1)
print("a pour forme normale de Smith H1=")
print(H1)
print("et matrice de transformation L1=")
print(L1)
print("et matrice de transformation U1=")
print(U1)
print("La matrice X2 = ")
print(X2)
print("a pour forme normale de Smith H2=")
print(H2)
print("et matrice de transformation L2=")
print(L2)
print("et matrice de transformation U2=")
print(U2)

print("\n$ Question 5")
print("La forme normale de Smith sur Q est ")
print(AQ)
print("La forme normale de Smith sur F2 est ")
print(AF2)
print("La forme normale de Smith sur F3 est ")
print(AF3)
print("La forme normale de Smith sur F5 est ")
print(AF5)

print("\n$ Question 6")
print("Votre fonction coincide avec celle de Sage ?")
print(test)


reset()
print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X = matrix(ZZ,[
      [ -2,  3,  3,  1],
      [  2, -1,  1, -3],
      [ -4,  0, -1, -4]])

# Code pour l'EXERCICE
def SageHNF(A):
    """
    Forme normale d'Hermite de SAGE avec la convention francaise :
    Les vecteurs sont ecrits en colonne.
    """
    m = A.nrows()
    n = A.ncols()
    Mm = identity_matrix(ZZ,m)[::-1,:]
    Mn = identity_matrix(ZZ,n)[::-1,:]
    AA = (Mm*A).transpose()
    HH, UU = AA.hermite_form(transformation=True)
    H = (HH*Mm).transpose()*Mn
    U = UU.transpose()*Mn
    assert(H-A*U==0)
    return H,U

H,U = SageHNF(X)
print("Forme normale d'Hermite")
print(H)
print("Matrice U")
print(U)

def indices_colonnes_nulles(H):
    n = H.ncols()
    indices = []
    for j in range(n):
        if H[:,j] == 0:
            indices.append(j)
    return indices

L =[] # liste des vecteurs d'une base
indices = indices_colonnes_nulles(H)
for i in indices:
    L.append(U[:,i])

# # Affichage des resultats

print("Le systeme a pour racine le module engendre par")
print(L)

reset()
print("""\
# ****************************************************************************
# IMAGE D'UNE MATRICE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A  = matrix(ZZ, [
           [ 15,  8, -9, 23,  -9],
           [ 22, 22,  7, -8,  20],
           [ 21, 18, -1, -7,  -3],
           [  3, -1,  0, 12, -16]])


# Code pour l'EXERCICE

def SageHNF(A):
    """
    Forme normale d'Hermite de SAGE avec la convention francaise :
    Les vecteurs sont ecrits en colonne.
    """
    m = A.nrows()
    n = A.ncols()
    Mm = identity_matrix(ZZ,m)[::-1,:]
    Mn = identity_matrix(ZZ,n)[::-1,:]
    AA = (Mm*A).transpose()
    HH, UU = AA.hermite_form(transformation=True)
    H = (HH*Mm).transpose()*Mn
    U = UU.transpose()*Mn
    assert(H-A*U==0)
    return H,U
#1
print("Mise sous forme normale d'Hermite")
H,U = SageHNF(A)
print(H)
print("Im(M) n'est pas égal à Z4, par exemple, on ne peut obtenir (1,1,0,1) avec la base de H")
print("On aurait x_1 = k, 2x_2 = 1, x_3 = 1, x_4 = 0 et x_5 = 1, or il n'existe pas d'entier tel que 2x_2 = 1")
print("""\
# 2.
""")
Z4 = ZZ^4
M = Z4.submodule(A.transpose())
test = (M == Z4)

# # Affichage des resultats

print("L'image de")
print(A)
print("est-elle egale a ZZ^4 ?")
print(test)

reset()
print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE NON-HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,[
           [ -6,  12,  6],
           [ 12, -16, -2],
           [ 12, -16, -2]])

b1 = vector(ZZ,[ -6, 4, 4])

PolF5.<x> = PolynomialRing(GF(5))

X2 = matrix(PolF5,[
           [ x + 1, 2,     4],
           [     1, x,     2],
           [     1, 1, x + 1]])

b2 = vector(PolF5,[ 3*x+2, 0, -1])

X3 = matrix(ZZ,[
    [2,-3,4,0],
    [4,4,53,-5]
])
b3 = vector(ZZ,[-4,2])

# Code pour l'EXERCICE


def solution_existe(a,b):
    for i in range(len(a)):
        if b[i]%a[i]!=0:
            return False
    
    for i in range(len(a)+1,len(b)):
        if b[i] != 0:
            return False
    return True     
       
def solve_nh(X,b):
    t = b[0].parent()
    S = X.smith_form()
    H = S[0]
    
    a = H.diagonal()
    
    for i in range(len(a)):
        if a[i]==0:
            a = a[:i]
            break
            
    L = S[1]
    C = S[2]
    n = C.nrows()
    b1 = L*b
    if not solution_existe(a,b1):
        return ("Pas de solution","Pas de solution")
        
    s_part = zero_vector(t,n)
    
    for i in range(len(a)):
        v = vector(C[:,i])
        q = b1[i]//a[i]
        s_part =s_part + v*q
        
    H = []
    for i in range(len(a),len(b)):
        H.append(C[:,i])
    return (s_part,H)

z1,H1 = solve_nh(X1,b1)

z2,H2 = solve_nh(X2,b2)

z3,H3 = solve_nh(X3,b3)

# # Affichage des resultats

print("Une solution particuliere de X1*z1 = b1 est")
print(z1)
print("les solutions du systeme homogene sont engendres par")
print(H1)
print("Une solution particuliere de X2*z2 = b2 est")
print(z2)
print("les solutions du systeme homogene sont engendrees par")
print(H2)
print("Une solution particuliere du systeme 3 est")
print(z3)
print("les solutions du systeme homogene sont engendres par")
print(H3)



reset()
print("""\
# ****************************************************************************
# STRUCTURE DU QUOTIENT
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])

# Code pour l'EXERCICE
S = A.smith_form()
H = S[0]
L = S[1]
C = S[2]
a = H.diagonal()
for i in range(len(a)):
    if a[i] == 0:
        a = a[:i]
        break
print("Inverse de L")
M = L.inverse()
print(M)
print("facteurs invariants")
print(a)
l = [vector(a[0]*M[:,0]),vector(a[1]*M[:,1]), vector(a[2]*M[:,2])]
print("l :")
print(l)
reponse = "(Z/15Z) + (Z/105Z) + (Z/630Z). Les vecteurs de l sont une base adaptée de N, on en déduit la structure suivante (Z/15Z)M1 + (Z/105Z)M2 + (Z/630Z)M3 où M1,M2,M3 sont colonnes de l'inverse de L, on en déduit que Z^3/N est isomorphe à (Z/15Z) + (Z/105Z) + (Z/630Z) "

# # Affichage des resultats

print("La structure de Z^3/N est")
print(reponse)

reset()
print("""\
# ****************************************************************************
# FACTEURS INVARIANTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])


# Code pour l'EXERCICE

rang = 0 # Par la structure de Z^3/N (cf. exo 84)
fact_inv = [15,105,630] # Idem, par l'exo 84
reponse = "Les exposants de M sont les multiples de 630, car M est isomorphe à Z/15Z + Z/105Z + Z/630Z. a*(1,1,1) = 0 si a est exposant, donc a doit être multiple de 15,105 et 630 étant donné que 15|105|630, on a facilement que les seuls exposants possibles sont les multiples de 630. Et réciproquement, les multiples de 630 vérifies 630k M = 0 "

# # Affichage des resultats

print("Le rang de Z^3 / N est")
print(rang)
print("Les facteurs invariants sont")
print(fact_inv)
print("Exposants ?")
print(reponse)

