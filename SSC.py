import time
import numpy as np
import matplotlib.pyplot as plt
import Synchrotron
import InverseCompton

delta = 14
dL = 1e27
R = 4.2e15
B = 0.2
z = 0

K = 1
n1 = 1.8
n2 = 2.7
gamma_min = 1
gamma_break = 1.1e5
gamma_max = 3.0e6
SpectrumType = 2
SpectrumPars = (K, n1, n2, gamma_min, gamma_break, gamma_max)

BList = np.logspace(-3, 3, 10)
KList = np.logspace(0, 5, 10)
nList = np.linspace(1.5, 4, 10)
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for K in KList]
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for n in nList]
#SpectrumParsList = [(K, n, gamma_min, gamma_max) for K in KList for n in nList]

nu = np.logspace(5, 25, 200)
#for B in BList:
#for SpectrumPars in SpectrumParsList:
tmin = time.time()
j_Synchrotron = Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
k_Synchrotron = Synchrotron.K(nu, B, SpectrumType, SpectrumPars)
#j_InverseCompton = InverseCompton.J(nu, B, R, SpectrumType, SpectrumPars, list(nu[j_Synchrotron > 1e-99]), list(j_Synchrotron[j_Synchrotron > 1e-99]), list(k_Synchrotron[j_Synchrotron > 1e-99]))
j_InverseCompton = InverseCompton.J(nu, B, R, SpectrumType, SpectrumPars, list(nu), list(j_Synchrotron), list(k_Synchrotron))
tmax = time.time()
print(f'SSC emission use {tmax-tmin}s')
plt.loglog(nu, nu*j_Synchrotron, 'r-', label=f'B={B}T')
plt.loglog(nu, nu*j_InverseCompton, 'b-', label=f'B={B}T')
plt.xlim(1e5, 1e25)
plt.ylim(1e-25, 1e-5)
plt.show()
