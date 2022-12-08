import numpy as np
import scipy.special as sc
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#### Solution exacte du probleme suivant
#  1/alpha \d_t T = 1/r \d_r (r \d_r T)    0 <= r < b
#  T(r,t) = f(t)  pour r = b
#  T(r,t) = 0     pour t = 0
# d'apres le livre Heat Conduction d'Ozisik (1980) p 204 ou 280, on trouve
#
# T(r,t) = 2 alpha/b * \sum_{n=1}^{\infty} beta_n J0(r beta_n)/J1(b beta_n) \int_0^t (f(tau) exp(-alpha*beta_n**2 (t-tau)) dtau) 
# ou (beta_n)_{1,\infty} sont les racines de J0(b x) = 0
#
# Formulation equivalente pour faire apparaitre le fait que T(b,t) = f(t)
#
# T(r,t) = f(t) - 2/b * \sum_{n=1}^{\infty} J0(r beta_n)/(beta_n*J1(b beta_n)) (f(0)exp(-alpha beta_n**2-t) + \int_0^t (f'(tau) exp(-alpha*beta_n**2 (t-tau)) dtau ))
#


#alpha = kappa/(rho*cp)
alpha = 0.001
# Rmax
b = 0.1
 
# Condition de dirichlet en r=b
f  = lambda t: 1500#t**3
fp = lambda t: 0 # 3*t**2

# Zero de la fonction J0(x*b)
# erreur ici comment tenir compte de b ??
M = 1000
beta_n = sc.jn_zeros(0, M)/b

# Calcul de la solution exacte
def sol(r,t):
    somme = 0.
    # for i in range(0,M):
    #     betan = beta_n[i]
    #     somme+= betan * sc.j0(r*betan)/sc.j1(b*betan) * integrate.quad(lambda tau: f(tau)*np.exp(-alpha*betan**2*(t-tau)), 0, t)[0]
    # sol = 2.*alpha/b * somme
    for i in range(0,M):
        betan = beta_n[i]
        somme+= sc.j0(r*betan)/(betan*sc.j1(b*betan)) * (f(0)*np.exp(-alpha*betan**2*t) + integrate.quad(lambda tau: fp(tau)*np.exp(-alpha*betan**2*(t-tau)), 0, t)[0])
    sol = f(t) - 2./b * somme

    return sol



# 100 linearly spaced numbers
r = np.linspace(0,b,100)
dt = 0.01
Tmax=100*dt

# fig = plt.figure() # initialise la figure
# line, = plt.plot([], []) 
# plt.xlim(0, b)
# plt.ylim(0, f(101*dt))

# def animate(i): 
#     t = i * dt
#     line.set_data(r, sol(r,t))
#     return line,
# ani = animation.FuncAnimation(fig, animate, frames=100, blit=True, interval=20, repeat=False)
# plt.show()

fsol = open("exactSol.dat", "w")
for i in xrange(0,100):
    fsol.write(str(r[i])+' '+str(sol(r,Tmax)[i])+"\n")
fsol.close()
