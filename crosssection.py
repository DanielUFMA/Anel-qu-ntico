#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.integrate import quad
import numpy as np
from scipy.special import gamma as gamma_function
from scipy.constants import e, c, alpha
import matplotlib.pyplot as plt

#constants

epsilon = 13.1          # dielectric constant of the material
gamma_C = 1           # donor impurity linewidth
nr = 3.2                # refractive index of semiconductor
flux = 0.0              # Phi in eqn 8 magnetic flux
R = 5.0e-9                # radius of the quantum ring in m
h = 1/2*np.pi            # Planck constant in eV
hbar = 1 #6.58e-16  #weird hack       # reduced Planck constant in eV
rho = np.linspace(0, 100,100)
nu = np.linspace(0, 100, 1000)
# Function that calculates the integrand
def func(rho):
    betai = gamma**2 / 2.0
    betaf = np.sqrt(1.0 + gamma**4 / 4.0)
    return (gamma * rho)**(betai + betaf) * np.exp(-0.5 * (gamma * rho)**2) * (gamma*rho**2)/2

def compute_matrix_element(gamma):
    betai = gamma**2 / 2.0
    betaf = np.sqrt(1.0 + gamma**4 / 4.0)
    integral = np.pi * quad(func, 0, np.infty)[0]
    first_sqrt  = (1.0 / (2.0**betai * gamma_function(betai + 1)))**0.5
    second_sqrt = (1.0 / (2.0**betaf * gamma_function(betaf + 1)))**0.5
    return R * gamma**2 / 2.0 / np.pi * first_sqrt * second_sqrt * integral

def compute_cross_section(hnu, gamma):
    # function that calculates the photoionisation cross section
    betai = gamma**2 / 2.0
    betaf = np.sqrt(1.0 + gamma**4 / 4.0)
    Ei = gamma**2 * (1.0 + betai) - gamma**4 / 2.0
    Ef = gamma**2 * (1.0 + betaf) - gamma**4 / 2.0
    delta = hbar*gamma_C / ((h*nu - Ef - Ei)**2 + (hbar * gamma_C)**2)
    matrix_element = compute_matrix_element(gamma)
    return nr / epsilon * 4.0 * np.pi / 3.0 * alpha * h*nu * abs(matrix_element)**2 * delta

#Plot


plt.figure()
for gamma in [1.0, 1.5, 2.0]:
    plt.plot(nu/h, (compute_cross_section(h*nu, gamma)))

plt.legend(['$\gamma = 1.0$', '$\gamma = 1.5$', '$\gamma = 2.0$'] )
plt.ylabel('$\\times$-section $\sigma$ (cm$^{2}$)')
plt.xlabel('Photon energy $h\\nu$ (meV)')
plt.show()
#plt.xlim(0,100)
