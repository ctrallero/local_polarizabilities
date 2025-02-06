import numpy as np
    
def V6(l1 = None,l2 = None,m = None): 
    V12 = (- 1) ** (m) * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * (16 / 231 * Wigner3j(l1,6,l2,0,0,0) * Wigner3j(l1,6,l2,- m,0,m) + 24 / 77 * Wigner3j(l1,4,l2,0,0,0) * Wigner3j(l1,4,l2,- m,0,m) + 10 / 21 * Wigner3j(l1,2,l2,0,0,0) * Wigner3j(l1,2,l2,- m,0,m)) + 1 / 7 * (l1 == l2)
    return V12
    
    return V12