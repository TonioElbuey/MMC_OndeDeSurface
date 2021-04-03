import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

# Données

nu = 0.25 # Coefficient de Poisson (ATTENTION INFLUENCE LES AUTRES)
c = 4.1e3 # Célérité ondes suivant X1
A = 10 # Amplitude (Arbitraire ?)
lamb = 5e3 # Longueur d'onde (Arbitraire ?)
k = 2*np.pi/lamb

H = 35e3 # Hauteur couche milieu a
p_a = 2.7 # Masse volumique milieu a
mu_a = 33e9 # Module de cisaillement milieu a
c_Ta = 3e3
c_La = 4.0e3

p_b = 3.3 # Masse volumique milieu b (ou début sujet)
mu_b = 67e9 # Module de cisaillement milieu b (ou début sujet)
c_Tb = 4.5e3
c_Lb = 5.2e3

b_a = k*np.sqrt((c**2/c_Ta**2) - 1)
b = k*np.sqrt(1 - (c**2/c_Tb**2))

N = 1000

# Solutions (partie réelle)

def u2_a(x1, x3, t):
    return ( A*np.cos(b_a*x3) + A*(mu_b*b/(mu_a*b_a))*np.sin(b_a*x3) )*np.cos(k*(x1 - c*t))

def u2_b(x1, x3, t):
    return A*np.exp(b*x3)*np.cos(k*(x1 - c*t))

# Simulation avec X3 variable

t = 0
x1 = 0
x3_a = np.linspace(H, 0, N)
x3_b = np.linspace(-H/3, 0, N)

onde_a = u2_a( x1, x3_a, t)
onde_b = u2_b( x1, x3_b, t)

plt.figure("Onde de Love", figsize = (6, 8))
plt.title("Onde de Love")
plt.xlabel("u2 (m)")
plt.ylabel("Profondeur X3 (m)")
plt.axhline(0, linestyle = '--', color = "k")
plt.axhline(H, linestyle = '--', color = "k")
plt.plot( onde_a,  x3_a, linestyle = '-', color = "r")
plt.plot( onde_b, x3_b, linestyle = '-', color = "b")
plt.tight_layout()
plt.savefig("Love_X3var.png")
plt.show()

# Simulation avec X1 variable et X3

fig = plt.figure("Onde de Love", figsize = (6, 8))
plt.title("Onde de Love")
plt.xlabel("u2 (m)")
plt.ylabel("X1 (m)")

t = 0
x1_min = 0
x1_max = 20e3
x1 = np.linspace(x1_min, x1_max, N)
onde_a, onde_b = [], []
line, = plt.plot([], [], color = 'r')
text = plt.gca().text(3., 100., '', fontsize = 10)
Nbis = 200

def init():
    global x1_min, x1_max
    plt.gca().set_xlim(-15, 15)
    plt.gca().set_ylim(x1_min, x1_max)

def update(i):

    global x1, line, onde_a, onde_b, t, text
    x3 = i # Choix profondeur
    x3_km = round(x3/1000, 1)

    if x3 < 0: # Dans la croute terrestre
        onde_b = u2_b( x1, x3, t)
        line.set_data(onde_b, x1)
    else: # Dans la couche superficielle
        onde_a = u2_a( x1, x3, t)
        line.set_data(onde_a, x1)

    text.set_text(f"Profondeur X3 = {x3_km} km")

anim = FuncAnimation(fig, update, np.linspace(-H/5, H, Nbis), init_func= init, interval = 500)
writer = PillowWriter(fps = 4)
anim.save("Love_X1X3var.gif", writer)
plt.show()