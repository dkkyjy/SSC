import time
import numpy as np
import matplotlib.pyplot as plt
import InverseCompton

delta = 10
dL = 1e27
R = 1e16
B = 1
z = 0

K = 1
n = 2
gamma_min = 1e2
gamma_max = 1e5
SpectrumType = 1
SpectrumPars = (K, n, gamma_min, gamma_max)

nu_Synchrotron = np.load('../Synchrotron/nu.npy')
nu = np.logspace(15, 25, 100)
for B in [1]:#, 10, 100]:
    j_Synchrotron = np.load(f'../Synchrotron/j_{B}.npy')
    k_Synchrotron = np.load(f'../Synchrotron/k_{B}.npy')

    tmin = time.time()
    j = InverseCompton.J(nu, B, R, SpectrumType, SpectrumPars, list(nu_Synchrotron[j_Synchrotron > 1e-99]), list(j_Synchrotron[j_Synchrotron > 1e-99]), list(k_Synchrotron[j_Synchrotron > 1e-99]))
    tmax = time.time()
    print(f'USING {tmax - tmin}s')
    plt.loglog(nu, nu*j, label=f'B={B}T')
plt.show()
