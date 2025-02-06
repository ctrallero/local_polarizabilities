import numpy as np
    
def V(l1 = None,l2 = None,m = None): 
    V12 = (- 1) ** m * 2 / 3 * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * Wigner3j(l1,2,l2,0,0,0) * Wigner3j(l1,2,l2,- m,0,m) + 1 / 3 * (l1 == l2)
    return V12
    
    return V12