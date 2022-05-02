import  numpy as np
import  matplotlib.pyplot as plt


a,b,c,d=1.5,0.05,0.48,0.005

#a=taux de reproduction intrinsèque des proies (constant, indépendant du nombre de prédateurs)
#b=taux de mortalité des proies dû aux prédateurs rencontrés
#c=taux de reproduction des prédateurs en fonction des proies rencontrées et mangées
#d=taux de mortalité intrinsèque des prédateurs (constant, indépendant du nombre de proies)

#x est le nombre de proies et y est le nombre de prédateurs (en fonction du temps) déterminés par les équations différentielles de Lotka Volterra
#ces équations sont approximées par la méthode d'Euler


#Variation de l'effectif des proies au cours du temps:
def proie(x,y):
    return a*x-b*x*y


#Variation de l'effectif des prédateurs au cours du temps:
def predateur(x,y):
    return -c*y+d*x*y


#Méthode d'Euler:
def Lotka_Volterra(x_0,y_0,tmin,tmax,h):
    liste_t, liste_x, liste_y=[tmin], [x_0], [y_0]
    t=tmin
    dy=y_0
    dx=x_0
    while t<=tmax:
        t+=h
        liste_t.append(t)
        dx+=proie(dx,dy)*h
        dy+=predateur(dx,dy)*h
        liste_x.append(dx)
        liste_y.append(dy)
    return liste_t,liste_x,liste_y


#Effectif des proies en fonction du temps:
def affichage_Lotka_Volterra_TX(x_0,y_0,tmin,tmax,h):
    T,X,Y=Lotka_Volterra(x_0,y_0,tmin,tmax,h)
    plt.plot(T,X)
    plt.title('Effectif des proies au cours du temps\nConditions initiales: 40 proies pour 10 prédateurs sur une durée de 50 ans')
    plt.xlabel('Temps (en année)')
    plt.ylabel('X = Effectif des proies')
    plt.show()

#affichage_Lotka_Volterra_TX(40,10,0,50,0.0005)


#Effectif des prédateurs en fonction du temps:
def affichage_Lotka_Volterra_TY(x_0,y_0,tmin,tmax,h):
    T,X,Y=Lotka_Volterra(x_0,y_0,tmin,tmax,h)
    plt.plot(T,Y)
    plt.title('Effectif des prédateurs au cours du temps\nConditions initiales: 40 proies pour 10 prédateurs sur une durée de 50 ans')
    plt.xlabel('T = Temps (en année)')
    plt.ylabel('Y = Effectif des prédateurs')
    plt.show()

#affichage_Lotka_Volterra_TY(40,10,0,50,0.0005)


#Effectifs des proies et des prédateurs en fonction du temps:
def affichage_Lotka_Volterra_TX_et_TY(x_0,y_0,tmin,tmax,h):
    T,X,Y=Lotka_Volterra(x_0,y_0,tmin,tmax,h)
    plt.plot(T,X)
    T,X,Y=Lotka_Volterra(x_0,y_0,tmin,tmax,h)
    plt.plot(T,Y)
    plt.title('Effectifs des prédateurs/proies au cours du temps\nConditions initiales: 40 proies pour 10 prédateurs sur une durée de 50 ans')
    plt.xlabel('T = Temps (en année)')
    plt.ylabel('Y = Effectif des prédateurs/proies')
    plt.show()

affichage_Lotka_Volterra_TX_et_TY(40,10,0,50,0.0005)


#Effectif des prédateurs en fonction de l'effectif des proies
def affichage_Lotka_Volterra_YX(x_0,y_0,tmin,tmax,h):
    T,X,Y=Lotka_Volterra(x_0,y_0,tmin,tmax,h)
    plt.plot(X,Y, label="X0=40 et Y0=10")
    T,X,Y=Lotka_Volterra(x_0+100,y_0+100,tmin,tmax,h)
    plt.plot(X,Y, label="X0=140 et Y0=110")
    plt.title('Effectif des proies en fonction de l-effectif des prédateurs\nConditions initiales : X0 proies pour Y0 prédateurs sur une durée de 50 ans et des relevés tous les 0.0005 ans')
    plt.xlabel('X = Effectif des proies')
    plt.ylabel('Y = Effectif des prédateurs')
    plt.legend(loc="upper right")
    plt.show()

affichage_Lotka_Volterra_YX(40,10,0,50,0.0005)