# Gráfico do potencial não-harmônico, axialmente simétrico com um núcleo centrífugo

import numpy as np
import matplotlib.pyplot as plt

R = 5  # nm
r = np.linspace(1.0e-15, 6 * R,100)
rho = r / R

def V(gamma):
    return gamma**4/4*(rho**2 + 1/rho**2) - gamma**4/2

for gamma in [1.0, 2.0, 5.0, 10.0]:
    sol = V(gamma)
    plt.plot(rho, sol, label = '$\gamma$ = ' + str(gamma))

plt.title("Potencial Anel quântico")
plt.legend()
plt.ylabel('$V(r)/E_{0}$')
plt.xlabel('r/R')
