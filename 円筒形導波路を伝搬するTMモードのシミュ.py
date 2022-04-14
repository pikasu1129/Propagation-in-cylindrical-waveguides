import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.axes_grid1
from scipy.special import jv

m = 0
n = 1

def djv(m,r, dr = 0.01):
    return (jv(m, r+dr) - jv(m, r))/dr

def ZERO_Bessel(m, r0 = 0, dr = 0.1, max_iter = 100, eps = 1e-10):
    r = r0
    while jv(m, r0)*jv(m, r) > 0:
        r = r + dr

    for i in range(max_iter):
        R = (r0 + r)/2
        if np.abs(R - r) < eps:
            break
        if jv(m, R) * jv(m, r) < 0:
            r0 = R
        if jv(m, R) * jv(m, r) > 0:
            r = R
        if i == max_iter - 1:
            print("doesn't calc")
    return R

def Xmn(m,n):
    Xmn = [0]
    for j in range(n):
        Xmn.append(ZERO_Bessel(m, r0 = Xmn[j] + 0.1))
    del Xmn[0]
    return Xmn

def TM(m, n, Xmn, theta, r, E = 1, epsilon = 1, mu = 1, omega = 1, t = 0, z = 0, flag = 1):
    kc = Xmn/a
    gamma = 1j
    Emn = 1j * E
    Ez = Emn * jv(m, kc*r)*np.exp(1j*(-1j*gamma*z + omega*t + flag*m*theta))
    Er = -gamma/kc**2*Emn*djv(m, kc*r)*np.exp(1j*(-1j*gamma*z + omega*t + flag*m*theta))
    Etheta = -flag*1j*gamma*m/(kc**2*r)*Emn*jv(m, kc*r)*np.exp(1j*(-1j*gamma*z + omega*t + flag*m*theta))
    Br = - (1j*mu*epsilon*omega/gamma)*Etheta
    Btheta = (1j*mu*epsilon*omega/gamma)*Er
    return Er.real, Etheta.real, Br.real, Btheta.real

def Cyliner_plot():
    fig = plt.figure(figsize = (6,4))
    axE = fig.add_subplot(121, polar = True)
    axE.streamplot(Theta, R, Etheta, Er, color = 'r', density = [1, 0.5])
    axE.set_ylim(0,a)
    axE.set_title('Electric field')
    axB = fig.add_subplot(122, polar = True)
    axB.streamplot(Theta, R, Btheta, Br, color = 'b', density = [1, 0.5])
    axB.set_ylim(0,a)
    axB.set_title('Magnetic field')
    plt.tight_layout()
    if m == 0:
        fig1 = plt.figure(figsize = (5,3))
        axE1 = fig1.add_subplot(121)
        axE1.plot(r, Er[:, 0], color = 'r', label = 'Er')
        axE1.plot(r, Etheta[:, 0], color = 'r', linestyle = 'dashed', label = r'$E_{\theta}$')
        axE1.set_ylabel('E')
        axE1.set_xlabel('r')
        axE1.set_xlim(0,a)
        axE1.legend()

        axB1 = fig1.add_subplot(122)
        axB1.plot(r, Br[:, 0], color = 'b', label = 'Br')
        axB1.plot(r, Btheta[:, 0], color = 'b', linestyle = 'dashed', label = r'$B_{\theta}$')
        axB1.set_ylabel('B')
        axB1.set_xlabel('r')
        axB1.set_xlim(0,a)
        axB1.legend()
        plt.tight_layout()

a = 1
theta = np.linspace(0, 2*np.pi, 360)
r = np.linspace(0.01, a, 100)
Theta, R = np.meshgrid(theta, r)

Xmn = Xmn(m,n)[n-1]
Er, Etheta, Br, Btheta = TM(m, n, Xmn, Theta, R, a)
Cyliner_plot()

plt.show()