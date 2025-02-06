# Wigner3j.m by David Terr, Raytheon, 6-17-04

# Calculates the Wigner 3j symbol using the Racah formula [1].

import numpy as np
from math import log, exp


def logfact(n):
    """
    Compute the logarithm of n!
    :param n: Input integer
    :return: log(n!)
    """
    a = 0
    for k in np.arange(1, n):
        a = a + log(n - k + 1)
    return a


def Wigner3j(j1, j2, j3, m1, m2, m3):
    """
    Function to compute the Wigner 3j symbol using the Racah formula
    :param j1, j2, j3: Angular momentum quantum numbers
    :param m1, m2, m3: Magnetic quantum numbers
    :return: Value of the Wigner 3j symbol
    """
    # Error checking
    if (2 * j1 != int(np.floor(2 * j1)) or 2 * j2 != int(np.floor(2 * j2)) or 2 * j3 != int(np.floor(2 * j3)) or
            2 * m1 != int(np.floor(2 * m1)) or 2 * m2 != int(np.floor(2 * m2)) or 2 * m3 != int(np.floor(2 * m3))):
        raise Exception('All parameters must be integers or half-integers.')
        return

    if j1 - m1 != int(np.floor(j1 - m1)):
        raise Exception('2*j1 and 2*m1 must have the same parity')
        return

    if j2 - m2 != int(np.floor(j2 - m2)):
        raise Exception('2*j2 and 2*m2 must have the same parity')
        return

    if j3 - m3 != int(np.floor(j3 - m3)):
        raise Exception('2*j3 and 2*m3 must have the same parity')
        return

    # Check the triangle condition: j3 must be between |j1 - j2| and j1 + j2
    if j3 > j1 + j2 or j3 < abs(j1 - j2):
        return 0.0

    # Check the range of magnetic quantum numbers
    if abs(m1) > j1:
        return 0.0

    if abs(m2) > j2:
        return 0.0

    if abs(m3) > j3:
        return 0.0

    # Compute the bounds for t
    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max(0, t1, t2)  # Minimum value of t
    tmax = min(t3, t4, t5)  # Maximum value of t

    wigner = 0.0  # Initialize the Wigner 3j symbol value

    # If m1 + m2 + m3 == 0, compute the Wigner 3j symbol
    if m1 + m2 + m3 == 0:
        for t in range(tmin, tmax + 1):
            # Compute the sign factor (-1)^t
            sign = (-1) ** t

            # Compute the logarithm of the denominator
            log_denominator = (
                    logfact(t) + logfact(t - t1) + logfact(t - t2) +
                    logfact(t3 - t) + logfact(t4 - t) + logfact(t5 - t)
            )

            # Compute the value of the current term
            term = sign * np.exp(-log_denominator)

            # Accumulate into wigner
            wigner += term

    sign = (-1.0) ** (j1 - j2 - m3)

    # Compute sqrt(numerator / denominator)
    log_sqrt = 0.5 * (logfact(j1 + j2 - j3) + logfact(j1 - j2 + j3) + logfact(-j1 + j2 + j3)
                      - logfact(j1 + j2 + j3 + 1) +
                      logfact(j1 + m1) + logfact(j1 - m1) +
                      logfact(j2 + m2) + logfact(j2 - m2) +
                      logfact(j3 + m3) + logfact(j3 - m3))

    # Compute the final value of the Wigner 3j symbol
    wigner = wigner * sign * np.exp(log_sqrt)

    return wigner

# Reference: Wigner 3j symbol entry in Eric Weinstein's Mathworld:
# http://mathworld.wolfram.com/Wigner3j-Symbol.html