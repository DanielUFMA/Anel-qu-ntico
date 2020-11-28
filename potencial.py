#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:27:31 2020

@author: daniel
"""
import numpy as np
import matplotlib.pyplot as plt
R = 5.0  # nm
r = np.linspace(10e-15, 10 * R, 100)
rho = r / R

def potencial(rho):
    V = gamma**4*(rho**2 + 1/rho**2) - (gamma**4)/2
    return V

Gamma = [(1.0, 'red'),
        (2.0, 'green'),
        (5.0, 'blue'),
        (10.0, 'black')]

fig, ax = plt.subplots()

for gamma in Gamma:
    ax.plot(rho, potencial(rho), color=gamma[1], label=r'$\gamma ='+str(gamma[0])+'$')

ax.legend()
ax.set_xlabel(r'$V(r)/E_{0}$')
ax.set_ylabel(r'$r/R$')