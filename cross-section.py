from scipy.integrste import quad
import numpy as np
from scipy.special import gamma as gamma_function
from scipy.constants import e, c
import matplolib.pyplot as plt

#function that calculates the integrand
def func(rho):
    betai = gamma**2/2.0
    betaf = np.sqrt(1.0 + gamma**4/4.0)
    return (gamma * rho)**(betai + betaf) * np.exp(-0.5 * gamma * rho**2))

def matrix_element(gamma):
