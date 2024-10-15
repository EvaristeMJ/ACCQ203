print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP12 : CRYPTANALYSE ET CRYPTOGRAPHIE A BASE DE RESEAUX                      #
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
# SAC A DOS (NE PAS TRAITER LA QUESTION 2)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

b = [356,278,417,27,132,464,521]
s = 1287

# Code pour l'EXERCICE

B = matrix.identity(len(b)+1,ZZ)
for i in range(len(b)):
    B[len(b),i] = b[i]
B[-1,-1] = -s

basis = (B.T).LLL()

x = []

for v in basis:
    if sum(v[i]*b[i] for i in range(len(b))) <= s and all(v[i]>= 0 for i in range(len(b))):
        x = v

# # Affichage des resultats

print("Le message est")
print(x)
print("Pour s'en assurer, on remarque que la somme obtenue avec ce vecteur est :", sum(x[i]*b[i] for i in range(len(b))), "On ne peut donc pas faire mieux.")

reset()
print("""\
# ****************************************************************************
# ATTAQUE DE WIENER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

N1 = 65946239999
e1 = 22022476093
N2 = 65946239999
e2 = 10865199773

phi_N1 = euler_phi(N1)
phi_N2 = euler_phi(N2)
E1 = mod(e1,phi_N1)
E2 = mod(e2,phi_N2)
s1 = E1^(-1)
s2 = E2^(-1)

B = floor(sqrt(sqrt(float(N1))))

# Code pour l'EXERCICE
ans  = (s1,"clé numéro 1")
if s1 < s2:
    ans = (s2, "clé numéro 2")
if s1 <= B and s2 <= B:
    ans = "Aucune ne devrait être utilisée car trop petite."



# # Affichage des resultats

print("Il vaut mieux utiliser la cle ",ans)

reset()
print("""\
# ****************************************************************************
# METHODE DE COPPERSMITH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

Pol.<x> = PolynomialRing(ZZ)

# Code pour l'EXERCICE

def Mat(f,m,B,N):
    
    d = f.degree()
    G = []
    for i in range(m+1):
        for j in range(i*d+1):
            temp = f^(m-i)
            g = x^j*N^i*temp
            g=g(B*x)
            G.append(g.padded_list(2*d+1))
    M = matrix(ZZ,G)
            
    return M

def Coppersmith(f,N):
    d = f.degree()
    m  = ceil(log(N,10)/d)
    B = ceil((N^(1/d))/(2*e))
    M = Mat(f,m,B,N)
    basis = (M).LLL()
    R = []
    i = 0
    while basis[i,:] == 0:
        i = i +1
    v = basis[i]
    h = sum((v[i]*(x)^i)/(B^i) for i in range(len(v)))
    
    return h.roots(ZZ)


# # Affichage des resultats

p=(x+1)*(x-2)*(x-3)*(x-29)
f = x^2-27*x-33
print(Coppersmith(p,10000))

reset()
print("""\
# ****************************************************************************
# MESSAGES STEREOTYPES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

bin=BinaryStrings()
N = 42564360034887861127
Pol.<x> = PolynomialRing(ZZ)
PolmodN.<y> = PolynomialRing(Integers(N))
e = 3
c = 12843085802751039909

# Code pour l'EXERCICE

mm=str(bin.encoding("08/06:00"))
mm = '0'*(8*ceil(len(mm)/8)-len(mm))+mm
mm = ZZ(mm,base=2)
f = PolmodN((mm+y)^e-c)
R = f.small_roots()
r = R[0]
m = mm+r
m = ZZ(m%N).str(2)
m = '0'*(8*ceil(len(m)/8)-len(m))+m
m = bin(m).decoding()

# # Affichage des resultats

print("Ce jour la, le message est")
print(m)

reset()
print("""\
# ****************************************************************************
# ALGORITHME DE BABAI (NE TRAITER QUE LA QUESTION 1)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

def Babai(B,t):
    # B est déjà LLL-réduite
    n = B.nrows()
    B_star,_ = B.gram_schmidt()
    b = t
    for j in range(n-1,-1,-1):
        u = b*B_star[j]/(B_star[j]*B_star[j])
        b = b- round(u)*B[j]
    return t-b

# # Affichage des resultats


reset()
print("""\
# ****************************************************************************
# CRYPTOSYSTEME GGH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

EE= matrix(ZZ,40, {(32, 32): 1, (22, 13): -2, (23, 8): -3,
       (15, 31): -1, (22, 37): -1,
       (13, 5): -4, (38, 20): 2, (4, 12): 3, (19, 22): -2, (15, 5): 2, (11,
       32): -1, (11, 10): 3, (1, 11): -4, (12, 33): 1, (0, 15): 1, (33, 17): 1,
       (7, 19): -1, (11, 1): -2, (7, 27): 3, (19, 32): -4, (22, 10): 2, (31,
       39): -4, (34, 9): 2, (36, 17): 2, (18, 17): 1, (14, 6): -2, (23, 14): 3,
       (23, 34): 2, (12, 11): -3, (0, 21): -3, (27, 22): -2, (4, 29): -3, (23,
       5): 1, (4, 6): -2, (24, 7): 2, (5, 38): -2, (33, 13): -1, (9, 35): 3,
       (18, 36): 1, (22, 5): 1, (24, 25): 3, (34, 31): 2, (6, 34): -3, (23,
       33): -4, (20, 37): -1, (38, 12): 2, (33, 0): -1, (4, 32): 3})
AA=10*identity_matrix(40)+EE
HH,U=AA.hermite_form(transformation = True)
cc = vector([-2, 0, 2, 0, 0, 1, -1, -1, -3, 0, 0, 2, -1, 13, 7, 2, 0, 2, 27, 2, 1,
       17, -2, 899, 50, 15, 11, 1098, 7, 2, -1, 10, -1, 2, 156, 15, 42, 8,
       525748584, 37])


# Code pour l'EXERCICE
bin = BinaryStrings()
def Babai(B,t):
    # B est déjà LLL-réduite
    B_star,_ = B.gram_schmidt()
    n = B_star.nrows()
    b = t
    v = vector([0 for i in range(len(t))])
    for j in range(n-1,-1,-1):
        u = b*B_star[j]/(B_star[j]*B_star[j])
        v += round(u)*B[j]
        b = b - (u-round(u))*B_star[j]-round(u)*B[j]
    return v

def StringToAscii(text):
    """
    Input : text String
    Output : list of ints
    """
    return [ZZ(ord(str(s))) for s in list(text)]

def AsciiToString(listascii):
    """
    Input : list of int
    Output : text String
    """
    return "".join([str(chr(a)) for a in listascii])

def StringToInts(text):
    """
    Input : text String
    Output : 
    """
    return [ZZ(str(t),base=2) for t in list(bin.encoding(text))]

def IntsToString(listints):
    return bin("".join([str(t) for t in listints])).decoding()

#print(m)
#print(bin(AsciiToString(m)).decoding())
#print(U*AA == HH)
mm = Babai(HH.LLL(),((AA.T)^(-1))*cc)

m = (U.T)^(-1)*mm

# # Affichage des resultats
print("Je n'ai pas réussi à obtenir le message en texte mais j'obtiens ce m après application de Babai : ",m)

