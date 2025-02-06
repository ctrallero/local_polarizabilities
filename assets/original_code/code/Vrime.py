import numpy as np
    
def Vprime(l1 = None,l2 = None,m = None): 
    Vprime = (- 1) ** m * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * (8 / 35 * Wigner3j(l1,4,l2,0,0,0) * Wigner3j(l1,4,l2,- m,0,m) + 4 / 7 * Wigner3j(l1,2,l2,0,0,0) * Wigner3j(l1,2,l2,- m,0,m) + 1 / 5 * Wigner3j(l1,0,l2,0,0,0) * Wigner3j(l1,0,l2,- m,0,m))
    return Vprime
    
    return Vprime