
# Wigner3j.m by David Terr, Raytheon, 6-17-04

# Compute the Wigner 3j symbol using the Racah formula [1].

import numpy as np
from math import *

def logfact(b):
    a=0
    for k in np.arange(1,b):
        a=a+log(b-k+1)
        
    return a
    
def Wigner3j(j1 = None,j2 = None,j3 = None,m1 = None,m2 = None,m3 = None): 
    # error checking
    if (2 * j1 != int(np.floor(2 * j1)) or 2 * j2 != int(np.floor(2 * j2)) or 2 * j3 != int(np.floor(2 * j3)) or 2 * m1 != int(np.floor(2 * m1)) or 2 * m2 != int(np.floor(2 * m2)) or 2 * m3 != int(np.floor(2 * m3))):
        raise Exception('All arguments must be integers or half-integers.')
        return
    
    if (j1 - m1 != int(np.floor(j1 - m1))):
        raise Exception('2*j1 and 2*m1 must have the same parity')
        return
    
    if (j2 - m2 != int(np.floor(j2 - m2))):
        raise Exception('2*j2 and 2*m2 must have the same parity')
        return
    
    if (j3 - m3 != int(np.floor(j3 - m3))):
        raise Exception('2*j3 and 2*m3 must have the same parity')
        return
    
    if j3 > j1 + j2 or j3 < np.abs(j1 - j2):
        wigner = 0
        return wigner
    
    if np.abs(m1) > j1:
        wigner = 0
        return wigner
    
    if np.abs(m2) > j2:
        wigner = 0
        return wigner
    
    if np.abs(m3) > j3:
        wigner = 0
        return wigner
    
    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2
    tmin = max(0,max(t1,t2))
    tmax = min(t3,min(t4,t5))
       
    wigner = 0
    if m1 + m2 + m3 == 0:
        if tmax > tmin:
            tarr = np.arange(tmin,tmax).reshape(-1)
            for t in tarr:
                wigner = wigner + np.exp(t * np.log(- 1) - (logfact(t) + logfact(t - t1) + logfact(t - t2) + logfact(t3 - t) + logfact(t4 - t) + logfact(t5 - t)))
                wigner = np.real(wigner)
        else:
            t=1
            wigner = wigner + np.exp(t * np.log(- 1) - (logfact(t) + logfact(t - t1) + logfact(t - t2) + logfact(t3 - t) + logfact(t4 - t) + logfact(t5 - t)))
            
    #      for t = tmin:tmax
#          wigner = wigner + (-1)^t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) ...
#               * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )
#      end
    
    wigner = np.real(np.exp(np.log(wigner) + (j1 - j2 - m3) * np.log(- 1) + 1 / 2 * (logfact(j1 + j2 - j3) + logfact(j1 - j2 + j3) + logfact(- j1 + j2 + j3) - logfact(j1 + j2 + j3 + 1) + logfact(j1 + m1) + logfact(j1 - m1) + logfact(j2 + m2) + logfact(j2 - m2) + logfact(j3 + m3) + logfact(j3 - m3))))
    return wigner
    
    #wigner = wigner * (-1)^(j1-j2-m3) ...
#    * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1)...
#       * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) );
    
    # Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Wigner3j-Symbol.html
