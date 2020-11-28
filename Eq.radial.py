#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:12:14 2020

@author: daniel
"""
import numpy as np
from scipy.special import genlaguerre, gamma
import matplotlib.pyplot as plt

m = 0 # magnetic number
flux = 0  # Phi in eqn 8
R = 5  # nm
r = np.linspace(0, 6 * R,100)
rho = r / R

def R0(n, gama):
    beta = np.sqrt((m - flux)**2 + gama**4/4)
    raiz = np.sqrt(gamma(n + 1) / ( 2**beta * gamma(n + beta + 1)))
    genlag =  genlaguerre(n, beta)(gama**2 * rho**2 / 2)
    asymbeh = (gama * rho)**beta *  np.exp(-gama**2 * rho**2 / 4)
    return 1/np.sqrt(2*np.pi) * gama/R * raiz * asymbeh * genlag

Gama = [(1.5, 'red'),
        (2.0, 'green'),
        (2.5, 'blue'),
        (3.0, 'black')]

fig, ax = plt.subplots()

for gama in Gama:
    sol = R0(0, gama[0])
    ax.plot(rho, sol, color=gama[1], label=r'$\gamma ='+str(gama[0])+'$')

ax.legend()
ax.set_xlabel(r'$R/r$')
ax.set_ylabel(r'$R_0(r)$')