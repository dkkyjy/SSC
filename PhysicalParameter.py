"""
    Follow the paper A&A 367, 809-825 (2001)
    The multifrequency emission of Mrk 501 from radio to TeV gamma-rays
"""

delta = 14
B = 0.2
R = 4.2e15 #cm
z = 
n2 = 
c = 3.00e10 #speed of light [cm / s]

def Nu2Gamma_break(nu_break):
    C = 3.7e6
    return np.sqrt((1+z)*nu_break / (C * B * delta))

def Nu2Gamma_peak(nu_peak):
    C = -3.25e5 * n2**3 + 1.67e6 * n2**2 - 3.62e6 * n2 + 4.65e6
    return np.sqrt((1+z)*nu_break / (C * B * delta))

def t2R(t_var):
    return c * t_var * delta / (1+z)
