print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP6 : BASES DE GROEBNER ET SYSTEMES POLYNOMIAUX MULTIVARIES                 #
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
#  FONCTIONS DE SAGEMATH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f = 2*x^2*y+7*z^3

# Code pour l'EXERCICE

print("x<y^2 est : " , x<y^2)
print("f.lt() renvoie ",f.lt())
print("f.lc() renvoie ",f.lc())
print("f.lm() envoie ", f.lm())

reponse  ="On commence par définir un ensemble de polynômes multivariés d'indéterminées x,y,z sur le corps des rationnels, muni de l'ordre lexicographique.\n Par exemple, x<y^2 renvoie False, car le degré de x^1 + y^0 + z^0 > x^0+y^2+z^0 avec l'ordre lexicographique. \n f.lt() renvoie le plus grand terme pour l'ordre dont on a muni nos polynômes (lexicographique ici). \n f.lc() renvoie le coefficient de ce terme dominant. \n f.lm() en renvoie le monôme dominant"
reponse = reponse + "\n Les différents ordres monomiaux de SageMath sont degrevlex, deglex, invlex,neglex, negdegrevlex, negdeglex, wdeglex,wdegrevlex,negwdegrevlex,degneglex."
reponse = reponse + "\n La présence d'un neg signifie qu'une condition < devient > pour la condition juste après negdeg => sur la condition de degré total, neglex => sur la comparaison lexicographique. La présence d'un lex signifie que c'est lexigographique. La présence d'un deg signifique qu'on fait d'abord une comparaison sur le degré total puis on passe au lexicographique en cas d'égalité. inv signifie qu'on fait marche l'ordre lex dans le sens inverse donc deg z, puis deg y puis deg x. w signifie qu'on ajoute des coefficients (weight) lors de la comparaison des degrés totaux (exemple : si on a un poids (1,2,3) alors deg(x^1 + y^3 + z^0)_w = 1 + 6 = 7 < deg(z^3)_w = 9 )"
# # Affichage des resultats

print("\n$1/ ", reponse)

reset()
print("""\
# ****************************************************************************
# DIVISION MULTIVARIEE
# ****************************************************************************
""")




# Donnees de l'enonce de l'exercice

MPol.<x,y> = PolynomialRing(QQ,2, order='lex')
f  = -x^7 + x^6*y + 2*x^5 - 2*x^4*y - 5*x^2 + 3*x*y^3 + 5*x*y + 11*y^3 + 10 
f1 = x*y^2+2*y^2
f2 = x^5+5

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    Q = [MPol(0)]*s
    p = f
    r = 0
    while p!=0:
        div = false
        for i in range(s):
            if p.lt()%F[i].lt() == 0:
                q = (p.lt())//(F[i].lt())
                Q[i] = Q[i]+q
                p = p - q*F[i]
                div = true
                break
        if not div:
            r = r + p.lt()
            p = p - p.lt()
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

# # Affichage des resultats

print("\n$ ",  myDivision(f,[f1,f2]))

reset()
print("""\
# ****************************************************************************
# BASE DE GROEBNER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x^2-y
f2 = x*y-z
f3 = z^4+x*y

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    Q = [MPol(0)]*s
    p = f
    r = 0
    if f.degree()==0:
        return Q,f
    while p!=0:
        div = false
        for i in range(s):
            if F[i] == 0:
                continue
            if p.lt()%F[i].lt() == 0:
                q = (p.lt())//(F[i].lt())
                Q[i] = Q[i]+q
                p = p - q*F[i]
                div = true
                break
        if not div:
            r = r + p.lt()
            p = p - p.lt()
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

def S(f,g):
    """
    Entrée : deux polynômes multivariés f,g dans K[x1,...,xn]
    Sortie : polynôme de syzygie de f et g
    """
    return (f.lt().lcm(g.lt())//f.lt())*f - (f.lt().lcm(g.lt())//g.lt())*g

def myGroebner(F):
    """
    Entrée : F liste de polynômes qui engendrent l'idéal J = <f1,...,fn>
    Sortie : Une base de Gröebner de J = <f1,...,fn>
    """
    G = F
    while True:
        L = []
        for i in range(len(G)):
            for j in range(len(G)):
                r = S(G[i],G[j]) 
                r = myDivision(r,G)[1]
                if r !=0:
                    L.append(r)
        G.extend(set(L)-set(G))
        if len(L) == 0:
            break
            
    return G

def estUneBaseGroebner(F,G):
    """
    Entrée : F une famille de polynômes, G la base de Groebner à teste
    Sortie : True si G est une base de Groebner de l'idéal engendré par F, False sinon
    """
    if Ideal(G) != Ideal(F):
        return False
    
    for i in range(len(G)):
        for j in range(i+1,len(G)):
            if myDivision(S(G[i],G[j]),G)[1] != 0:
                return False
    
    return True

def reduction(G):
    """
    Entrée : Base de Groebner
    Sortie : Base de Groebner après une étape dans le processus de réduction
    """
    GG = G.copy()
    i = 0
    while i < len(GG):
        g = GG[i]
        F = GG.copy()
        F.remove(g)
        _,r = myDivision(g,F)
        if r!=0:
            GG[i] = r*(r.lc())^(-1)
            i = i + 1
        else:
            GG.remove(g)
    if set(GG) == set(G):
        return G
    return reduction(GG)
        

def myRedGroebner(F):
    G = myGroebner(F)
    Gred = reduction(G)
    return Gred

G = myGroebner([f1,f2,f3])
Gred = myRedGroebner([f1,f2,f3])
I = Ideal([f1,f2,f3])
# # Affichage des resultats

print("\n$1/ ",G)
print("Est-ce une base de Groebner ?",estUneBaseGroebner([f1,f2,f3],G))
print("\n$2/ ",Gred)
print("Est-ce une base Groebner minimale et réduite ?",set(I.groebner_basis())==set(Gred))



reset()
print("""\
# ****************************************************************************
# APPARTENANCE A UN IDEAL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x*y-y^2
f2 = x^3-z^2
I = Ideal([f1,f2])
f = -4*x^2*y^2*z^2 + y^6 + 3*z^5

# Code pour l'EXERCICE
def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    Q = [MPol(0)]*s
    p = f
    r = 0
    if f.degree()==0:
        return Q,f
    while p!=0:
        div = false
        for i in range(s):
            if p.lt()%F[i].lt() == 0:
                q = (p.lt())//(F[i].lt())
                Q[i] = Q[i]+q
                p = p - q*F[i]
                div = true
                break
        if not div:
            r = r + p.lt()
            p = p - p.lt()
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r


test1 = f in I
test2 = myDivision(f,[f1,f2])[1] == 0


reset()
print("""\
# ****************************************************************************
# RESOLUTION D'UN SYSTEME
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


MPol.<x,y> = PolynomialRing(QQ,2,order="lex") # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = (y^2+6)*(x-1) - y*(x^2 + 1)
g = (x^2+6)*(y-1) - x*(y^2 + 1)
 

# Code pour l'EXERCICE
J = Ideal([f,g])
base = J.groebner_basis() # Vous pouvez utiliser la fonction adhoc de sage
          # pour calculer la base Groebner
h = base[2].univariate_polynomial() # On peut remarquer que l'un des polynômes est en uniquement en y (le polynôme 3)
racines_y =[racine_y for racine_y,_ in h.roots()]
racines  = []

# On distringue les cas selon les racines de h le polynôme en y
for pol in base:
    for racine_y in racines_y:
        sub_poly = pol.subs(y = racine_y).univariate_polynomial() # On évalue en les racines, les polynômes sont alors en x
        if sub_poly == 0:
            continue
        for racine_x,_ in sub_poly.roots():
            if (racine_x,racine_y) not in racines:
                racines.append((racine_x,racine_y)) # On calcule alors les tuples de zéros (x,y)


Gf = implicit_plot(f,(x,0,6),(y,0,6),color='red') 
Gg = implicit_plot(g,(x,0,6),(y,0,6),color='blue')  
Gp = point2d(racines,color='green')

# # Affichage des resultats

print("\n$1/  Une base de Groebner de [f,g] est", base)
print("\n$2/  Les valeurs de y sont", racines_y)
print("\n$3/  Les valeurs de (x,y) sont", racines)
print("\n$4/")
show(Gf+Gg+Gp)

reset()
print("""\
# ****************************************************************************
# OPTIMISATION
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


MPol.<x,y,lamb> = PolynomialRing(QQ,3,order="invlex") # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = x^2*y  - 2*x*y + y + 1
g = x^2 + y^2 - 1


# Code pour l'EXERCICE
df = f.gradient()
dg = g.gradient()
syst = [df[i]-lamb*dg[i] for i in range(2)] + [x^2+y^2-1]
J = Ideal(syst)
base = J.groebner_basis()
hx = base[-1]
racines_x = [racine_x for racine_x,_ in hx.univariate_polynomial().roots(RR)]
racines_x_y = [(racine_x,racine_y) for racine_x in racines_x for racine_y,_ in base[1].subs(x = racine_x).univariate_polynomial().roots(RR)]
racines = [(racine_x,racine_y) for (racine_x,racine_y) in racines_x_y for l,_  in base[0].subs(x=racine_x,y=racine_y).univariate_polynomial().roots(RR)]

# # Affichage des resultats

print("\n$1/  On doit resoudre le systeme", syst)
print("\n$2/  dont une base de Groebner est", base)
print("Ici, on utilise une base de Groebner obtenue en utilisant un ordre invlex, on remarque que cela permet de directement trouver des racines de x")
print("\n$4/  Les valeurs de (x,y) sont",racines)
var('t')
p = f.subs(x=cos(t),y=sin(t))
plot(p(t),(t,0,10))



reset()
print("""\
# ****************************************************************************
# MANIPULATIONS ALGEBRIQUES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice
MPol.<c,s,u,v> = PolynomialRing(QQ,4,order="degrevlex")
# On pose c=cos theta s = sin theta, u est le u(theta) de l'énoncé, idem pour v
f = c+s - u 
g = c^2-s^2 + 2*c*s - v
h = c^2+s^2-1 
J = Ideal([f,g,h])
# Code pour l'EXERCICE

r = (s^6).reduce(J)

# # Affichage des resultats

print("\n$1/ L'expression de sin(theta)^6 en fonction de u et v est ")
print(r)

reset()
print("""\
# ****************************************************************************
# OVALES DE DESCARTES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,OM,PM> = PolynomialRing(QQ,4,order="invlex")


# Code pour l'EXERCICE

f = x^2+y^2-OM^2
g= (x-1)^2+y^2-PM^2
h = OM+2*PM-3

J = Ideal([f,g,h])
eq = J.groebner_basis()[-1] # On remarque que le dernier élément est une équation cartésienne

var('x,y')
eq = eq(x,y,0,0)
Geq = implicit_plot(eq,(x,-3,5),(y,-5,5),color='red') 

# # Affichage des resultats

print("\n$ L'équation est ",eq," = 0")
show(Geq)