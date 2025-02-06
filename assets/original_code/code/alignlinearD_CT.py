import numpy as np
import time as time
from scipy import interpolate
from scipy.integrate import solve_ivp

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from sympy.physics.wigner import wigner_3j
from sympy import *
init_printing(use_latex=False)

from alignlinearmultipulses import *

from matplotlib import pyplot

def alignlinearD(molecule_dictionary): 
    #alignlinearD(filename,tau,I, tau2, I2)
    # inputs are file with molecular configurations, pulse durations tau and
    # tau2 in fs and intensities of the two beams in TW/cm^2
    
    mol_name = molecule_dictionary['molecule_name']
    rotcon = molecule_dictionary['B']
    centrifugal = molecule_dictionary['centrifugal']
    
    T= molecule_dictionary['T']
    Jmax = molecule_dictionary['Jmax']
    reneven = molecule_dictionary['reneven']
    renodd = molecule_dictionary['renodd']
    torp1 = molecule_dictionary['tau']
    I01 = molecule_dictionary['I'] * 10 ** 12
    Tfree1 = molecule_dictionary['Tfree1']
    Tfree2 = molecule_dictionary['Tfree2']
    Tfree3 = molecule_dictionary['Tfree3']
    torp2 = molecule_dictionary['tau2']
    I02 = molecule_dictionary['I02'] * 10 ** 12
    torp3 = molecule_dictionary['tau3']
    I03 = molecule_dictionary['I03'] * 10 ** 12
    NOPlaser = molecule_dictionary['NOPlaser']
    stepfree = molecule_dictionary['stepfree']    
    alphaper = molecule_dictionary['alphaper']
    alphapar = molecule_dictionary['alphapar']
    
    
    fname1 = mol_name + ' T=' + str(T) + 'K Jmax=' + str(Jmax) + ' tau1=' + \
            str(torp1) + ' I01=' + '%.2g'%I01 + ' tau2=' + str(torp2) +    \
            ' I02=' + '%.2g'%I02 + ' tau3=' + str(torp2) + ' I03=' +       \
            '%.2g'%I02 + ' .dat'    

    ## change to a.u
    lambda_ = 0.8
    
    hbar = 6.6 * 10 ** - 34 / 2 / np.pi
    me = 9.1 * 10 ** - 31
    epsion = 8.85 * 10 ** - 12
    qe = 1.6e-19
    es = qe / np.sqrt(4 * np.pi * epsion)
    t0 = hbar ** 3 / me / es ** 4
    w = 2 * np.pi * 3 / lambda_ * 10 ** 14
    w0 = 2 * np.pi * 3 / lambda_ * 10 ** 14 * t0
    c = 300000000.0
    a0 = 5.3e-11
    E01 = np.sqrt(I01) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2
    E02 = np.sqrt(I02) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2
    E03 = np.sqrt(I03) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2
    KB = 1.38 * 10 ** (- 23) * T
    kb = KB / (4.3597 * 10 ** - 18)
    pulsed1 = torp1 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)
    pulsedfree1 = Tfree1 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)
    pulsed2 = torp2 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)
    pulsedfree2 = Tfree2 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)
    pulsed3 = torp3 * 10 ** - 15 / (2.41888 * 10 ** - 17)
    pulsedfree3 = Tfree3 * 10 ** - 15 / (2.41888 * 10 ** - 17)
    stepfree = stepfree * 10 ** - 15 / (2.418884326505 * 10 ** - 17)
    #Tsep=Tsep*10^-15/(2.41888*10^-17);
    
    pulsed = np.zeros(3)
    pulsed[1-1] = pulsed1
    pulsed[2-1] = pulsed2
    pulsed[3-1] = pulsed3
    ## Electric field
    NOPfree1 = int(np.floor(pulsedfree1 / stepfree)) # Number of points for free evolution after pulse 1
    NOPfree2 = int(np.floor(pulsedfree2 / stepfree)) # Number of points for free evolution after pulse 2
    NOPfree3 = int(np.floor(pulsedfree3 / stepfree)) # Number of points for free evolution after pulse 3
    NOPfree = np.zeros(3, dtype=int)
    NOPfree[1-1] = NOPfree1
    NOPfree[2-1] = NOPfree2
    NOPfree[3-1] = NOPfree3
    t = np.zeros((3,NOPlaser))
    V00 = np.zeros((3,NOPlaser))
    V02 = np.zeros((3,NOPlaser))
    a1 = 2 * np.log(2) / (pulsed1 ** 2)
    t1 = np.linspace(0,2 * pulsed1,NOPlaser)
    
    t[1-1,:] = t1
    E1 = np.multiply(E01,np.exp(- a1 * (t1 - pulsed1) ** 2))
    a2 = 2 * np.log(2) / (pulsed2 ** 2)
    t2 = np.linspace(0,2 * pulsed2,NOPlaser)
    
    t[2-1,:] = t2
    E2 = np.multiply(E02,np.exp(- a2 * (t2 - pulsed2) ** 2))
    a3 = 2 * np.log(2) / (pulsed3 ** 2)
    t3 = np.linspace(0,2 * pulsed3,NOPlaser)
    
    t[3-1,:] = t3
    E3 = np.multiply(E03,np.exp(- a3 * (t3 - pulsed3) ** 2))
    V00[1-1,:] = - 0.25 * alphaper * E1 ** 2
    V02[1-1,:] = - 0.25 * E1 ** 2 * (alphapar - alphaper)
    V00[2-1,:] = - 0.25 * alphaper * E2 ** 2
    V02[2-1,:] = - 0.25 * E2 ** 2 * (alphapar - alphaper)
    V00[3-1,:] = - 0.25 * alphaper * E3 ** 2
    V02[3-1,:] = - 0.25 * E3 ** 2 * (alphapar - alphaper)
    # alphaper=alphaper/(0.53^3);
    
    rotcon = rotcon / (8065.54429 * 27.2114)
    centrifugal = centrifugal / (8065.54429 * 27.2114)
    
    # Determine the initial number of thermal J states
    abc = 0
    patitionf = 0
    sumw = 0
    while (np.exp(- ((abc + 1) * abc * rotcon - centrifugal * abc ** 2 * (abc + 1) ** 2) / kb) * (2 * abc + 1)) > 1e-05:

        abc = abc + 1 # This is the initial number of J states

    # Calculation of the partition function
    # mac is the total numer of J levels
    for x in np.arange(1,abc + 1).reshape(-1):
        if (np.mod(x - 1,2) == 0):
            macj = reneven
        else:
            macj = renodd
        patitionf = patitionf + macj * np.exp(- ((x - 1) * x * rotcon - centrifugal * (x - 1) ** 2 * (x) ** 2) / kb) * (2 * x - 1)
    
    for mac in np.arange(1,100).reshape(-1):
        if (np.mod(mac - 1,2) == 0):
            vinodj = reneven
        else:
            vinodj = renodd
        if sumw < 0.999:
            # determine the number of initial J states 'mac' until the sum of the weight of different J states converges to bigger than 0.999.
            sumw = sumw + vinodj * np.exp(- ((mac - 1) * mac * rotcon - centrifugal * (mac - 1) ** 2 * (mac) ** 2) / kb) * (2 * (mac - 1) + 1) / patitionf
        else:
            break
    
    
    J = np.arange(0,mac+1,1)
    # Ej=J.*(J+1)*rotcon;
    Ej = np.multiply(J,(J + 1)) * rotcon - centrifugal * J ** 2.0 * (J + 1) ** 2
    tfree1 = np.linspace(2 * pulsed1,pulsedfree1,NOPfree1)
    tfree2 = np.linspace(2 * pulsed2,pulsedfree2,NOPfree2)
    tfree3 = np.linspace(2 * pulsed3,pulsedfree3,NOPfree3)
    tfree = np.zeros((3,np.max(NOPfree)))
    tfree[1-1,np.arange(1-1,NOPfree[1-1])] = tfree1
    tfree[2-1,np.arange(1-1,NOPfree[2-1])] = tfree2
    tfree[3-1,np.arange(1-1,NOPfree[3-1])] = tfree3
    #rotational energy levels
# Totalstates=mac
    cba = 0
    fullsize = NOPlaser * 3 + NOPfree1 + NOPfree2 + NOPfree3
    cos8thetaaverage = 1j*np.zeros(fullsize)
    cos6thetaaverage = 1j*np.zeros(fullsize)
    cos4thetaaverage = 1j*np.zeros(fullsize)
    cos2thetaaverage = 1j*np.zeros(fullsize)
    k = 0
    for x in np.arange(1,mac + 2).reshape(-1):
        for y in np.arange(1,x+1).reshape(-1):
            k = k + 1
    
    jsize = k
    D2 = 1j*np.zeros((fullsize,jsize))
    D4 = 1j*np.zeros((fullsize,jsize))
    D6 = 1j*np.zeros((fullsize,jsize))
    D8 = 1j*np.zeros((fullsize,jsize))
    symm = np.zeros(jsize)
    Eout = np.zeros(jsize)
    ii = 0
    jsize
    for a in np.arange(1,mac+2).reshape(-1): # a is the current J level, from matlab a=1:mac+1
        print(str(a) + '/' + str(mac))
        if (np.mod(a - 1,2) == 0):
            renj = reneven
        else:
            renj = renodd
        m = np.arange(0,a,1) # from matlab it was m=0:1:a-1
        # python needs to have more than one element in an array to loop. 
        #It will not execute the last one
        earr = np.arange(1,a+1) # from matlab it was m=0:1:a-1
        if renj != 0:
            for e in earr:
                ii = ii + 1
                print('mJ state ' + str(m[e-1]))
                if (m[e-1] == 0):
                    a2,a4,a6,a8 = alignlinearmultipulses(Jmax,a - 1,0,rotcon,centrifugal,V00,V02,t,tfree,NOPlaser,NOPfree,fullsize,pulsed)
                    cos2thetaaverage = cos2thetaaverage + renj * np.exp(- Ej[a-1] / kb) / patitionf * a2
                    cos4thetaaverage = cos4thetaaverage + renj * np.exp(- Ej[a-1] / kb) / patitionf * a4
                    cos6thetaaverage = cos6thetaaverage + renj * np.exp(- Ej[a-1] / kb) / patitionf * a6
                    cos8thetaaverage = cos8thetaaverage + renj * np.exp(- Ej[a-1] / kb) / patitionf * a8
                    D2[:,ii-1] = 1 / 2 * (3 * a2 - 1)
                    D4[:,ii-1] = 1 / 8 * (35 * a4 - 30 * a2 + 3)
                    D6[:,ii-1] = 1 / 16 * (231 * a6 - 315 * a4 + 105 * a2 - 5)
                    D8[:,ii-1] = 1 / 128 * (6435 * a8 - 12012 * a6 + 6930 * a4 - 1260 * a2 + 35)
                    symm[ii-1] = renj
                    Eout[ii-1] = Ej[a-1]
                else:
                    b2,b4,b6,b8 = alignlinearmultipulses(Jmax,a - 1,m[e-1],rotcon,centrifugal,V00,V02,t,tfree,NOPlaser,NOPfree,fullsize,pulsed)
                    cos2thetaaverage = cos2thetaaverage + 2 * renj * np.exp(- Ej[a-1] / kb) / patitionf * b2
                    cos4thetaaverage = cos4thetaaverage + 2 * renj * np.exp(- Ej[a-1] / kb) / patitionf * b4
                    cos6thetaaverage = cos6thetaaverage + 2 * renj * np.exp(- Ej[a-1] / kb) / patitionf * b6
                    cos8thetaaverage = cos8thetaaverage + 2 * renj * np.exp(- Ej[a-1] / kb) / patitionf * b8
                    D2[:,ii-1] = 1 / 2 * (3 * b2 - 1)
                    D4[:,ii-1] = 1 / 8 * (35 * b4 - 30 * b2 + 3)
                    D6[:,ii-1] = 1 / 16 * (231 * b6 - 315 * b4 + 105 * b2 - 5)
                    symm[ii-1] = 2 * renj
                    Eout[ii-1] = Ej[a-1]
            cba = cba + 1
            #strcat('J=',num2str(a-1),'...','...', 'E=',num2str(Ej(a)),'...',num2str(cba),'/',num2str(mac+1))
    
    ## plot
    tfreeout1 = tfree1 * 2.41888 * 10 ** - 5
    
    tout1 = t1 * 2.41888 * 10 ** - 5
    tfreeout2 = (tfree2 + pulsedfree1) * 2.41888 * 10 ** - 5
    
    tout2 = (t2 + pulsedfree1) * 2.41888 * 10 ** - 5
    tfreeout3 = (tfree3 + pulsedfree1 + pulsedfree2) * 2.41888 * 10 ** - 5
    
    tout3 = (t3 + pulsedfree1 + pulsedfree2) * 2.41888 * 10 ** - 5
    time = np.concatenate((tout1,tfreeout1,tout2,tfreeout2,tout3,tfreeout3))
    
    cos2thetaaverage = np.transpose(cos2thetaaverage.real)
    cos4thetaaverage = np.transpose(cos4thetaaverage.real)
    cos6thetaaverage = np.transpose(cos6thetaaverage.real)
    cos8thetaaverage = np.transpose(cos8thetaaverage.real)
    # sin2thetaaverage=1-cos2thetaaverage;
    
    sin2thetacos2theta = cos2thetaaverage - cos4thetaaverage
    sin2thetacos4theta = cos4thetaaverage - cos6thetaaverage
    sin2thetacos6theta = cos6thetaaverage - cos8thetaaverage
    # cos2thetaaverage=[real(cos2thetaaveragein2)',real(cos2thetaaverage2)'];
    
    return cos2thetaaverage,D2,D4,D6,D8,time,Eout,symm


molecule_dictionary = {
    'molecule_name': 'N2', 
    'T':5,                # Temperature K
    'Jmax': 5,            # Max J value
    'B':1.9895,             # B cm-1
    'centrifugal':0.0,                # Centrifugal distorsion 
    'alphaper':9.8,         # Bohr^3
    'alphapar':15,          # Bohr^3
    'reneven':0.6667,         # g even J
    'renodd':0.3333,          # g off J
    'tau': 50,              # first pulse duration in fs
    'I': 3,                 # peak intensity in TW/cm2
    'tau2': 50,             # pulse duration in fs
    'I02': 3,                # pulse 2 peak intensity in TW/cm2
    'Tfree1': 9000,              # evolution after first pulse  in fs
    'Tfree2': 10,              # evolution after second pulse  in fs
    'tau3': 10,                # pulse duration of third pulse in fs
    'Tfree3': 10,              # evolution after third pulse in fs
    'I03': 0.0,                 # Intesity of third pulse (W/cm2)
    'NOPlaser': 25,            # field points in fs
    'stepfree': 20             # step size in fs
    }

# =============================================================================
mol_name = molecule_dictionary['molecule_name']
T = molecule_dictionary['T']
Jmax = molecule_dictionary['Jmax']
# reneven = molecule_dictionary['reneven']
# renodd = molecule_dictionary['renodd']
torp1 = molecule_dictionary['tau']
I01 = molecule_dictionary['I'] * 10 ** 12
# Tfree1 = molecule_dictionary['Tfree1']
# Tfree2 = molecule_dictionary['Tfree2']
# Tfree3 = molecule_dictionary['Tfree3']
torp2 = molecule_dictionary['tau2']
I02 = molecule_dictionary['I02'] * 10 ** 12
torp3 = molecule_dictionary['tau3']
I03 = molecule_dictionary['I03'] * 10 ** 12
# NOPlaser = molecule_dictionary['NOPlaser']
# stepfree = molecule_dictionary['stepfree']    
# alphaper = molecule_dictionary['alphaper']
# alphapar = molecule_dictionary['alphapar']
# 
# =============================================================================
cos2thetaaverage,D2,D4,D6,D8,time,Eout,symm = alignlinearD(molecule_dictionary)



fname1 = mol_name + ' T=' + str(T) + 'K Jmax=' + str(Jmax) + ' tau1=' + \
        str(torp1) + ' I01=' + '%.2g'%I01 + ' tau2=' + str(torp2) +    \
        ' I02=' + '%.2g'%I02 + ' tau3=' + str(torp2) + ' I03=' +       \
        '%.2g'%I02 + ' .dat'