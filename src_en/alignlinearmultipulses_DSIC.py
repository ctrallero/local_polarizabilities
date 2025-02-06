import numpy as np
from scipy.integrate import solve_ivp

from sympy import *
from Wigner3j import Wigner3j

init_printing(use_latex=False)


# Calculate the matrix element <Yl1m|cos^2 theta|Yl2m>
def V(l1=None, l2=None, m=None):
    V12 = (-1) ** m * 2 / 3 * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * \
          N(Wigner3j(l1, 2, l2, 0, 0, 0)) * N(Wigner3j(l1, 2, l2, -m, 0, m)) + \
          1 / 3 * (l1 == l2)
    return V12


# Calculate the matrix element <Yl1m|cos^6 theta|Yl2m>
def V6(l1=None, l2=None, m=None):
    V12 = (-1) ** m * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * \
          (16 / 231 * N(Wigner3j(l1, 6, l2, 0, 0, 0)) * N(Wigner3j(l1, 6, l2, -m, 0, m)) + \
           24 / 77 * N(Wigner3j(l1, 4, l2, 0, 0, 0)) * N(Wigner3j(l1, 4, l2, -m, 0, m)) + \
           10 / 21 * N(Wigner3j(l1, 2, l2, 0, 0, 0)) * N(Wigner3j(l1, 2, l2, -m, 0, m))) + \
          1 / 7 * (l1 == l2)

    return V12


# Calculate the matrix element <Yl1m|cos^8 theta|Yl2m>
def V8(l1=None, l2=None, m=None):
    V12 = (-1) ** m * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * \
          (128 / 6435 * N(Wigner3j(l1, 8, l2, 0, 0, 0)) * N(Wigner3j(l1, 8, l2, -m, 0, m)) + \
           64 / 495 * N(Wigner3j(l1, 6, l2, 0, 0, 0)) * N(Wigner3j(l1, 6, l2, -m, 0, m)) + \
           48 / 143 * N(Wigner3j(l1, 4, l2, 0, 0, 0)) * N(Wigner3j(l1, 4, l2, -m, 0, m)) + \
           40 / 99 * N(Wigner3j(l1, 2, l2, 0, 0, 0)) * N(Wigner3j(l1, 2, l2, -m, 0, m))) + \
          1 / 9 * (l1 == l2)
    return V12


# Calculate the matrix element <Yl1m|cos^4 theta|Yl2m>
def Vprime(l1, l2, m):
    Vprime = (-1) ** m * np.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * \
             (8 / 35 * N(Wigner3j(l1, 4, l2, 0, 0, 0)) * N(Wigner3j(l1, 4, l2, -m, 0, m)) + \
              4 / 7 * N(Wigner3j(l1, 2, l2, 0, 0, 0)) * N(Wigner3j(l1, 2, l2, -m, 0, m)) + \
              1 / 5 * N(Wigner3j(l1, 0, l2, 0, 0, 0)) * N(Wigner3j(l1, 0, l2, -m, 0, m)))
    return Vprime


# # Define the field function for the ODE solver
# def fieldode(tfield, coeff, t, n, V00, V02, Ej, Vcos2theta):
#     dcoeff = np.zeros(n)
#     V00 = np.interp(tfield, t, V00)  # Interpolate V00
#     V02 = np.interp(tfield, t, V02)  # Interpolate V02
#     E = np.zeros((n, n))
#
#     for jj in np.arange(1, n + 1):
#         E[jj - 1, jj - 1] = Ej[jj - 1]  # Construct the energy matrix
#
#     return -1j * (np.matmul(E + V00 * np.eye(n) + V02 * Vcos2theta, coeff))

################################################
################################################
# Modified code: In addition to H_rot and H_ind, H_self is also included
def fieldode(tfield, coeff, t, n, V00, V02, V04, Ej, Vcos2theta, Vcos4theta):
    V00 = np.interp(tfield, t, V00)  # Interpolate V00
    V02 = np.interp(tfield, t, V02)  # Interpolate V02
    V04 = np.interp(tfield, t, V04)  # Interpolate V04
    E = np.zeros((n, n))

    for jj in np.arange(1, n + 1):
        E[jj - 1, jj - 1] = Ej[jj - 1]  # Construct the energy matrix

    return -1j * (np.matmul(E + V00 * np.eye(n) + V02 * Vcos2theta + V04 * Vcos4theta, coeff))


################################################
################################################

# Calculate molecular alignment under multiple pulses
def alignlinearmultipulses(Jmax=None, J0=None, m0=None, rotcon=None, centrifugal=None, V00=None, V02=None, V04=None,
                           t=None, tfree=None, NOPlaser=None, NOPfree=None, fullsize=None, pulsed=None):
    # Determine the size of the interaction matrix
    # If the initial J0 is even (odd), Jmax should also be even (odd)
    if (np.mod(J0, 2) == 0):
        if (np.mod(Jmax, 2) != 0):
            Jmax = Jmax - 1
    else:
        if (np.mod(Jmax, 2) == 0):
            Jmax = Jmax - 1

    # Remove transitions that may cause m > J and remove unnecessary transitions from the matrix elements
    n = int(np.floor((Jmax - m0) / 2 + 1))
    J = np.arange(Jmax, 0 - 1, -2)  # 0-1 because arange does not include the last index
    J = np.flipud(J[0:n])

    # Create the <Yl1m|cos^2 theta|Yl2m> matrix: Selection rules: deltaJ=0,2,-2
    Vcos2theta = np.zeros((n, n))
    for j in np.arange(1, n + 1).reshape(-1):
        Vcos2theta[j - 1, j - 1] = V(J[j - 1], J[j - 1], m0)
        if j + 1 <= n:
            Vcos2theta[j + 1 - 1, j - 1] = V(J[j + 1 - 1], J[j - 1], m0)

    Vcos2theta = Vcos2theta + np.transpose(np.tril(Vcos2theta, -1))

    # Create the <Yl1m|cos^4 theta|Yl2m> matrix: Selection rules: deltaJ=0,2,-2,4,-4
    Vcos4theta = np.zeros((n, n))
    for l in np.arange(1, n + 1).reshape(-1):
        Vcos4theta[l - 1, l - 1] = Vprime(J[l - 1], J[l - 1], m0)
        if l + 1 <= n - 1:
            Vcos4theta[l + 1 - 1, l - 1] = Vprime(J[l + 1 - 1], J[l - 1], m0)
        if l + 2 <= n - 1:
            Vcos4theta[l + 2 - 1, l - 1] = Vprime(J[l + 2 - 1], J[l - 1], m0)

    Vcos4theta = Vcos4theta + np.transpose(np.tril(Vcos4theta, -1))

    # Create the <Yl1m|cos^6 theta|Yl2m> matrix: Selection rules: deltaJ=0,2,-2,4,-4,6,-6
    Vcos6theta = np.zeros((n, n))
    for l in np.arange(1, n + 1).reshape(-1):
        Vcos6theta[l - 1, l - 1] = V6(J[l - 1] - 1, J[l - 1] - 1, m0)
        if l + 1 <= n - 1:
            Vcos6theta[l + 1 - 1, l - 1] = V6(J[l + 1 - 1], J[l - 1], m0)
        if l + 2 <= n - 1:
            Vcos6theta[l + 2 - 1, l - 1] = V6(J[l + 2 - 1], J[l - 1], m0)
        if l + 3 <= n - 1:
            Vcos6theta[l + 3 - 1, l - 1] = V6(J[l + 3 - 1], J[l - 1], m0)

    Vcos6theta = Vcos6theta + np.transpose(np.tril(Vcos6theta, -1))

    # Create the <Yl1m|cos^8 theta|Yl2m> matrix: Selection rules: deltaJ=0,2,-2,4,-4,6,-6,8,-8
    Vcos8theta = np.zeros((n, n))
    for l in np.arange(1, n + 1).reshape(-1):
        Vcos8theta[l - 1, l - 1] = V8(J[l - 1], J[l - 1], m0)
        if l + 1 <= n - 1:
            Vcos8theta[l + 1 - 1, l - 1] = V8(J[l + 1 - 1], J[l - 1], m0)
        if l + 2 <= n - 1:
            Vcos8theta[l + 2 - 1, l - 1] = V8(J[l + 2 - 1], J[l - 1], m0)
        if l + 3 <= n - 1:
            Vcos8theta[l + 3 - 1, l - 1] = V8(J[l + 3 - 1], J[l - 1], m0)
        if l + 4 <= n - 1:
            Vcos8theta[l + 4 - 1, l - 1] = V8(J[l + 3 - 1], J[l - 1], m0)

    Vcos8theta = Vcos8theta + np.transpose(np.tril(Vcos8theta, -1))

    # Calculate rotational energy
    Ej = np.multiply(J, (J + 1)) * rotcon - centrifugal * J ** 2.0 * (J + 1) ** 2

    coeffinit = np.transpose((J == J0))
    cos2theta = 1j * np.zeros(fullsize)
    cos4theta = 1j * np.zeros(fullsize)
    cos6theta = 1j * np.zeros(fullsize)
    cos8theta = 1j * np.zeros(fullsize)

    numpulse = 3  # Number of pulses
    for orz in np.arange(0, numpulse).reshape(-1):
        # sol = solve_ivp(fieldode, [0, 2 * pulsed[orz]], y0=np.complex_(coeffinit), dense_output=True,
        #                 args=(t[orz, :], n, V00[orz, :], V02[orz, :], Ej, Vcos2theta))

        ################################################
        ################################################
        sol = solve_ivp(fieldode, [0, 2 * pulsed[orz]], y0=np.complex_(coeffinit), dense_output=True,
                        args=(
                            t[orz, :], n, V00[orz, :], V02[orz, :], V04[orz, :], Ej, Vcos2theta, Vcos4theta,
                        ))
        ################################################
        ################################################

        tsol = np.linspace(0, 2 * pulsed[orz], 50)
        Y = sol.sol(tsol)
        T = tsol
        Y = Y.transpose()  # To match MATLAB ODE
        m1 = T.size
        cos2thetaprime = 1j * np.zeros(m1)
        cos4thetaprime = 1j * np.zeros(m1)
        cos6thetaprime = 1j * np.zeros(m1)
        cos8thetaprime = 1j * np.zeros(m1)
        for s in np.arange(1, m1 + 1).reshape(-1):
            cos2thetaprime[s - 1] = np.dot(Y[s - 1, :], np.matmul(Vcos2theta, np.conj(Y[s - 1, :])))
            cos4thetaprime[s - 1] = np.dot(Y[s - 1, :], np.matmul(Vcos4theta, np.conj(Y[s - 1, :])))
            cos6thetaprime[s - 1] = np.dot(Y[s - 1, :], np.matmul(Vcos6theta, np.conj(Y[s - 1, :])))
            cos8thetaprime[s - 1] = np.dot(Y[s - 1, :], np.matmul(Vcos8theta, np.conj(Y[s - 1, :])))

        pulsed_time_array_ind_left = (orz) * NOPlaser + (orz != 0) * NOPfree[0] + (orz == 2) * NOPfree[1]
        pulsed_time_array_ind_right = pulsed_time_array_ind_left + NOPlaser
        cos2theta[pulsed_time_array_ind_left:pulsed_time_array_ind_right] = np.interp(t[orz, :], T, cos2thetaprime)
        cos4theta[pulsed_time_array_ind_left:pulsed_time_array_ind_right] = np.interp(t[orz, :], T, cos4thetaprime)
        cos6theta[pulsed_time_array_ind_left:pulsed_time_array_ind_right] = np.interp(t[orz, :], T, cos6thetaprime)
        cos8theta[pulsed_time_array_ind_left:pulsed_time_array_ind_right] = np.interp(t[orz, :], T, cos8thetaprime)
        for mn in np.arange(0, NOPfree[orz]).reshape(-1):
            C = 1j * np.zeros(n)
            for j in np.arange(0, n).reshape(-1):
                C[j] = Y[m1 - 1, j] * np.exp(-1j * Ej[j] * (tfree[orz, mn] - 2 * pulsed[orz]))
            pulsed_time_array_ind = (orz + 1) * NOPlaser + (orz != 0) * NOPfree[0] + (orz == 2) * NOPfree[1]
            cos4theta[pulsed_time_array_ind + mn] = np.dot(np.conj(C), np.dot(Vcos4theta, C))
            cos2theta[pulsed_time_array_ind + mn] = np.dot(np.conj(C), np.dot(Vcos2theta, C))
            cos6theta[pulsed_time_array_ind + mn] = np.dot(np.conj(C), np.dot(Vcos6theta, C))
            cos8theta[pulsed_time_array_ind + mn] = np.dot(np.conj(C), np.dot(Vcos8theta, C))
        coeffinit = C

    return cos2theta, cos4theta, cos6theta, cos8theta