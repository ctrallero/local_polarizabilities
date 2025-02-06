import numpy as np
    
def V8(l1 = None,l2 = None,m = None): 
    V12 = (- 1) ** (m) * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * (128 / 6435 * Wigner3j(l1,8,l2,0,0,0) * Wigner3j(l1,8,l2,- m,0,m) + 64 / 495 * Wigner3j(l1,6,l2,0,0,0) * Wigner3j(l1,6,l2,- m,0,m) + 48 / 143 * Wigner3j(l1,4,l2,0,0,0) * Wigner3j(l1,4,l2,- m,0,m) + 40 / 99 * Wigner3j(l1,2,l2,0,0,0) * Wigner3j(l1,2,l2,- m,0,m)) + 1 / 9 * (l1 == l2)
    return V12
    
    return V12