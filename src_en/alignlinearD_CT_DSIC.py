# -*- coding: utf-8 -*-
# @Time: 2024/12/17 20:46
# @Software: PyCharm

import matplotlib.pyplot as plt

import warnings

from alignlinearD_CT import alignlinearD

warnings.simplefilter(action='ignore', category=FutureWarning)

from sympy import *

init_printing(use_latex=False)

from alignlinearmultipulses_DSIC import *


def alignlinearD_DSIC(molecule_dictionary):
    # Function description: Calculate the alignment behavior of molecules under laser pulses
    # Input parameters: Molecule configuration file, pulse width tau and tau2 (unit: fs), laser intensity (unit: TW/cm^2)

    mol_name = molecule_dictionary['molecule_name']  # Molecule name
    rotcon = molecule_dictionary['B']  # Rotational constant (unit: cm^-1)
    centrifugal = molecule_dictionary['centrifugal']  # Centrifugal distortion constant

    T = molecule_dictionary['T']  # Temperature (unit: K)
    Jmax = molecule_dictionary['Jmax']  # Maximum rotational quantum number J
    reneven = molecule_dictionary['reneven']  # Weight for even J states
    renodd = molecule_dictionary['renodd']  # Weight for odd J states
    torp1 = molecule_dictionary['tau']  # First pulse width (unit: fs)
    I01 = molecule_dictionary['I'] * 10 ** 12  # Peak intensity of the first pulse (unit: W/cm^2)
    Tfree1 = molecule_dictionary['Tfree1']  # Free evolution time after the first pulse (unit: fs)
    Tfree2 = molecule_dictionary['Tfree2']  # Free evolution time after the second pulse (unit: fs)
    Tfree3 = molecule_dictionary['Tfree3']  # Free evolution time after the third pulse (unit: fs)
    torp2 = molecule_dictionary['tau2']  # Second pulse width (unit: fs)
    I02 = molecule_dictionary['I02'] * 10 ** 12  # Peak intensity of the second pulse (unit: W/cm^2)
    torp3 = molecule_dictionary['tau3']  # Third pulse width (unit: fs)
    I03 = molecule_dictionary['I03'] * 10 ** 12  # Peak intensity of the third pulse (unit: W/cm^2)
    NOPlaser = molecule_dictionary['NOPlaser']  # Number of time points in the laser field
    stepfree = molecule_dictionary['stepfree']  # Free evolution step size (unit: fs)
    alphaper = molecule_dictionary['alphaper']  # Perpendicular polarizability (unit: Bohr^3)
    alphapar = molecule_dictionary['alphapar']  # Parallel polarizability (unit: Bohr^3)

    ## Convert to atomic units (a.u.)
    lambda_ = 0.8  # Laser wavelength (unit: μm)

    hbar = 6.6 * 10 ** - 34 / 2 / np.pi  # Reduced Planck constant (unit: J·s)
    me = 9.1 * 10 ** - 31  # Electron mass (unit: kg)
    epsion = 8.85 * 10 ** - 12  # Vacuum permittivity (unit: F/m)
    qe = 1.6e-19  # Electron charge (unit: C)
    es = qe / np.sqrt(4 * np.pi * epsion)  # Electrostatic energy of the electron
    t0 = hbar ** 3 / me / es ** 4  # Time unit conversion factor
    w = 2 * np.pi * 3 / lambda_ * 10 ** 14  # Laser frequency (unit: Hz)
    # w0 = 2 * np.pi * 3 / lambda_ * 10 ** 14 * t0  # Frequency converted to atomic units
    c = 300000000.0  # Speed of light (unit: m/s)
    a0 = 5.3e-11  # Bohr radius (unit: m)
    E01 = np.sqrt(I01) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2  # Electric field strength of the first pulse (atomic units)
    E02 = np.sqrt(I02) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2  # Electric field strength of the second pulse (atomic units)
    E03 = np.sqrt(I03) * 8.5 * 10 ** - 10 * lambda_ * me * c / qe * w * a0 / 27.2  # Electric field strength of the third pulse (atomic units)
    KB = 1.38 * 10 ** (- 23) * T  # Boltzmann constant (unit: J/K)
    kb = KB / (4.3597 * 10 ** - 18)  # Boltzmann constant converted to atomic units
    pulsed1 = torp1 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)  # First pulse width converted to atomic units
    pulsedfree1 = Tfree1 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)  # Free evolution time after the first pulse converted to atomic units
    pulsed2 = torp2 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)  # Second pulse width converted to atomic units
    pulsedfree2 = Tfree2 * 10 ** - 15 / (2.418884326505 * 10 ** - 17)  # Free evolution time after the second pulse converted to atomic units
    pulsed3 = torp3 * 10 ** - 15 / (2.41888 * 10 ** - 17)  # Third pulse width converted to atomic units
    pulsedfree3 = Tfree3 * 10 ** - 15 / (2.41888 * 10 ** - 17)  # Free evolution time after the third pulse converted to atomic units
    stepfree = stepfree * 10 ** - 15 / (2.418884326505 * 10 ** - 17)  # Free evolution step size converted to atomic units
    # Tsep = Tsep * 10^-15 / (2.41888 * 10^-17);  # Commented code, possibly for other time unit conversions

    ################################################
    ################################################
    lambda_dsic = 0.05  # Dipole self-interaction coupling constant
    ################################################
    ################################################

    pulsed = np.zeros(3)  # Array for pulse widths
    pulsed[1 - 1] = pulsed1  # First pulse width
    pulsed[2 - 1] = pulsed2  # Second pulse width
    pulsed[3 - 1] = pulsed3  # Third pulse width

    ## Electric field setup
    NOPfree1 = int(np.floor(pulsedfree1 / stepfree))  # Number of time points for free evolution after the first pulse
    NOPfree2 = int(np.floor(pulsedfree2 / stepfree))  # Number of time points for free evolution after the second pulse
    NOPfree3 = int(np.floor(pulsedfree3 / stepfree))  # Number of time points for free evolution after the third pulse
    NOPfree = np.zeros(3, dtype=int)  # Array for free evolution time points
    NOPfree[1 - 1] = NOPfree1  # Number of time points for free evolution after the first pulse
    NOPfree[2 - 1] = NOPfree2  # Number of time points for free evolution after the second pulse
    NOPfree[3 - 1] = NOPfree3  # Number of time points for free evolution after the third pulse

    t = np.zeros((3, NOPlaser))  # Time array
    V00 = np.zeros((3, NOPlaser))  # Induced potential V00 array
    V02 = np.zeros((3, NOPlaser))  # Induced potential V02 array

    # Electric field for the first pulse
    a1 = 2 * np.log(2) / (pulsed1 ** 2)  # Width parameter for the Gaussian pulse
    t1 = np.linspace(0, 2 * pulsed1, NOPlaser)  # Time points for the first pulse
    t[1 - 1, :] = t1  # Time points for the first pulse
    E1 = np.multiply(E01, np.exp(- a1 * (t1 - pulsed1) ** 2))  # Electric field strength for the first pulse

    # Electric field for the second pulse
    a2 = 2 * np.log(2) / (pulsed2 ** 2)  # Width parameter for the Gaussian pulse
    t2 = np.linspace(0, 2 * pulsed2, NOPlaser)  # Time points for the second pulse
    t[2 - 1, :] = t2  # Time points for the second pulse
    E2 = np.multiply(E02, np.exp(- a2 * (t2 - pulsed2) ** 2))  # Electric field strength for the second pulse

    # Electric field for the third pulse
    a3 = 2 * np.log(2) / (pulsed3 ** 2)  # Width parameter for the Gaussian pulse
    t3 = np.linspace(0, 2 * pulsed3, NOPlaser)  # Time points for the third pulse
    t[3 - 1, :] = t3  # Time points for the third pulse
    E3 = np.multiply(E03, np.exp(- a3 * (t3 - pulsed3) ** 2))  # Electric field strength for the third pulse

    ################################################
    ################################################
    V04 = np.zeros((3, NOPlaser))  # Induced potential V04 array

    # Calculate induced potentials V00, V02, V04
    V00[1 - 1, :] = - 0.25 * alphaper * (E1 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * (alphaper ** 2) * (E1 ** 2))
    V02[1 - 1, :] = - 0.25 * (alphapar - alphaper) * (E1 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * 2 * alphaper * (alphapar - alphaper) * (E1 ** 2))
    V04[1 - 1, :] = 0.25 * (lambda_dsic ** 2) * ((alphapar - alphaper) ** 2) * (E1 ** 2)

    V00[2 - 1, :] = - 0.25 * alphaper * (E2 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * (alphaper ** 2) * (E2 ** 2))
    V02[2 - 1, :] = - 0.25 * (alphapar - alphaper) * (E2 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * 2 * alphaper * (alphapar - alphaper) * (E2 ** 2))
    V04[2 - 1, :] = 0.25 * (lambda_dsic ** 2) * ((alphapar - alphaper) ** 2) * (E2 ** 2)

    V00[3 - 1, :] = - 0.25 * alphaper * (E3 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * (alphaper ** 2) * (E3 ** 2))
    V02[3 - 1, :] = - 0.25 * (alphapar - alphaper) * (E3 ** 2) + \
                    (0.25 * (lambda_dsic ** 2) * 2 * alphaper * (alphapar - alphaper) * (E3 ** 2))
    V04[3 - 1, :] = 0.25 * (lambda_dsic ** 2) * ((alphapar - alphaper) ** 2) * (E3 ** 2)
    ################################################
    ################################################

    # Convert rotational constant and centrifugal distortion constant to atomic units
    rotcon = rotcon / (8065.54429 * 27.2114)
    centrifugal = centrifugal / (8065.54429 * 27.2114)

    # Determine the number of initial thermal J states
    abc = 0
    patitionf = 0  # Partition function
    sumw = 0  # Sum of weights
    while (np.exp(- ((abc + 1) * abc * rotcon - centrifugal * abc ** 2 * (abc + 1) ** 2) / kb) * (2 * abc + 1)) > 1e-05:
        abc = abc + 1  # Number of initial J states

    # Calculate the partition function
    for x in np.arange(1, abc + 1).reshape(-1):
        if np.mod(x - 1, 2) == 0:
            macj = reneven
        else:
            macj = renodd
        patitionf = patitionf + macj * np.exp(- ((x - 1) * x * rotcon - centrifugal * (x - 1) ** 2 * (x) ** 2) / kb) * (
                2 * x - 1)

    # Determine the number of initial J states until the sum of weights converges to 0.999
    for mac in np.arange(1, 100).reshape(-1):
        if np.mod(mac - 1, 2) == 0:
            vinodj = reneven
        else:
            vinodj = renodd
        if sumw < 0.999:
            sumw = sumw + vinodj * np.exp(
                - ((mac - 1) * mac * rotcon - centrifugal * (mac - 1) ** 2 * (mac) ** 2) / kb) * (
                           2 * (mac - 1) + 1) / patitionf
        else:
            break

    J = np.arange(0, mac + 1, 1)  # Array of J states
    Ej = np.multiply(J, (J + 1)) * rotcon - centrifugal * J ** 2.0 * (J + 1) ** 2  # Rotational energy

    # Free evolution time points
    tfree1 = np.linspace(2 * pulsed1, pulsedfree1, NOPfree1)
    tfree2 = np.linspace(2 * pulsed2, pulsedfree2, NOPfree2)
    tfree3 = np.linspace(2 * pulsed3, pulsedfree3, NOPfree3)
    tfree = np.zeros((3, np.max(NOPfree)))
    tfree[1 - 1, np.arange(1 - 1, NOPfree[1 - 1])] = tfree1
    tfree[2 - 1, np.arange(1 - 1, NOPfree[2 - 1])] = tfree2
    tfree[3 - 1, np.arange(1 - 1, NOPfree[3 - 1])] = tfree3

    # Initialize variables
    cba = 0
    fullsize = NOPlaser * 3 + NOPfree1 + NOPfree2 + NOPfree3  # Total number of time points
    cos8thetaaverage = 1j * np.zeros(fullsize)  # <cos^8 theta> average
    cos6thetaaverage = 1j * np.zeros(fullsize)  # <cos^6 theta> average
    cos4thetaaverage = 1j * np.zeros(fullsize)  # <cos^4 theta> average
    cos2thetaaverage = 1j * np.zeros(fullsize)  # <cos^2 theta> average
    k = 0
    for x in np.arange(1, mac + 2).reshape(-1):
        for y in np.arange(1, x + 1).reshape(-1):
            k = k + 1

    jsize = k  # Number of J states
    D2 = 1j * np.zeros((fullsize, jsize))  # <cos^2 theta> matrix
    D4 = 1j * np.zeros((fullsize, jsize))  # <cos^4 theta> matrix
    D6 = 1j * np.zeros((fullsize, jsize))  # <cos^6 theta> matrix
    D8 = 1j * np.zeros((fullsize, jsize))  # <cos^8 theta> matrix
    symm = np.zeros(jsize)  # Symmetry weights
    Eout = np.zeros(jsize)  # Energy output
    ii = 0

    # Iterate over J states
    for a in np.arange(1, mac + 2).reshape(-1):
        print(str(a) + '/' + str(mac))
        if (np.mod(a - 1, 2) == 0):
            renj = reneven
        else:
            renj = renodd
        m = np.arange(0, a, 1)
        earr = np.arange(1, a + 1)
        if renj != 0:
            for e in earr:
                ii = ii + 1
                print('mJ state ' + str(m[e - 1]))
                if (m[e - 1] == 0):
                    # print(Jmax, a - 1, 0, rotcon, centrifugal, V00, V02, t, tfree, NOPlaser, NOPfree, fullsize, pulsed)
                    a2, a4, a6, a8 = alignlinearmultipulses(Jmax, a - 1, 0, rotcon, centrifugal, V00, V02, V04, t,
                                                            tfree,
                                                            NOPlaser, NOPfree, fullsize, pulsed)
                    cos2thetaaverage = cos2thetaaverage + renj * np.exp(- Ej[a - 1] / kb) / patitionf * a2
                    cos4thetaaverage = cos4thetaaverage + renj * np.exp(- Ej[a - 1] / kb) / patitionf * a4
                    cos6thetaaverage = cos6thetaaverage + renj * np.exp(- Ej[a - 1] / kb) / patitionf * a6
                    cos8thetaaverage = cos8thetaaverage + renj * np.exp(- Ej[a - 1] / kb) / patitionf * a8
                    D2[:, ii - 1] = 1 / 2 * (3 * a2 - 1)
                    D4[:, ii - 1] = 1 / 8 * (35 * a4 - 30 * a2 + 3)
                    D6[:, ii - 1] = 1 / 16 * (231 * a6 - 315 * a4 + 105 * a2 - 5)
                    D8[:, ii - 1] = 1 / 128 * (6435 * a8 - 12012 * a6 + 6930 * a4 - 1260 * a2 + 35)
                    symm[ii - 1] = renj
                    Eout[ii - 1] = Ej[a - 1]
                else:
                    b2, b4, b6, b8 = alignlinearmultipulses(Jmax, a - 1, m[e - 1], rotcon, centrifugal, V00, V02, V04,
                                                            t,
                                                            tfree, NOPlaser, NOPfree, fullsize, pulsed)
                    cos2thetaaverage = cos2thetaaverage + 2 * renj * np.exp(- Ej[a - 1] / kb) / patitionf * b2
                    cos4thetaaverage = cos4thetaaverage + 2 * renj * np.exp(- Ej[a - 1] / kb) / patitionf * b4
                    cos6thetaaverage = cos6thetaaverage + 2 * renj * np.exp(- Ej[a - 1] / kb) / patitionf * b6
                    cos8thetaaverage = cos8thetaaverage + 2 * renj * np.exp(- Ej[a - 1] / kb) / patitionf * b8
                    D2[:, ii - 1] = 1 / 2 * (3 * b2 - 1)
                    D4[:, ii - 1] = 1 / 8 * (35 * b4 - 30 * b2 + 3)
                    D6[:, ii - 1] = 1 / 16 * (231 * b6 - 315 * b4 + 105 * a2 - 5)
                    symm[ii - 1] = 2 * renj
                    Eout[ii - 1] = Ej[a - 1]
            cba = cba + 1

    ## Plotting
    tfreeout1 = tfree1 * 2.41888 * 10 ** - 5
    tout1 = t1 * 2.41888 * 10 ** - 5
    tfreeout2 = (tfree2 + pulsedfree1) * 2.41888 * 10 ** - 5
    tout2 = (t2 + pulsedfree1) * 2.41888 * 10 ** - 5
    tfreeout3 = (tfree3 + pulsedfree1 + pulsedfree2) * 2.41888 * 10 ** - 5
    tout3 = (t3 + pulsedfree1 + pulsedfree2) * 2.41888 * 10 ** - 5
    time = np.concatenate((tout1, tfreeout1, tout2, tfreeout2, tout3, tfreeout3))

    cos2thetaaverage = np.transpose(cos2thetaaverage.real)

    return cos2thetaaverage, D2, D4, D6, D8, time, Eout, symm, mol_name


# # Molecule parameter dictionary in the source code
# molecule_dictionary = {
#     'molecule_name': 'N2',
#     'T': 5,  # Temperature (unit: K)
#     'Jmax': 5,  # Maximum J value
#     'B': 1.9895,  # Rotational constant (unit: cm^-1)
#     'centrifugal': 0.0,  # Centrifugal distortion constant
#     'alphaper': 9.8,  # Perpendicular polarizability (unit: Bohr^3)
#     'alphapar': 15,  # Parallel polarizability (unit: Bohr^3)
#     'reneven': 0.6667,  # Weight for even J states
#     'renodd': 0.3333,  # Weight for odd J states
#     'tau': 50,  # First pulse width (unit: fs)
#     'I': 3,  # Peak intensity of the first pulse (unit: TW/cm^2)
#     'tau2': 50,  # Second pulse width (unit: fs)
#     'I02': 3,  # Peak intensity of the second pulse (unit: TW/cm^2)
#     'Tfree1': 9000,  # Free evolution time after the first pulse (unit: fs)
#     'Tfree2': 10,  # Free evolution time after the second pulse (unit: fs)
#     'tau3': 10,  # Third pulse width (unit: fs)
#     'Tfree3': 10,  # Free evolution time after the third pulse (unit: fs)
#     'I03': 0.0,  # Peak intensity of the third pulse (unit: W/cm^2)
#     'NOPlaser': 25,  # Number of time points in the laser field
#     'stepfree': 20  # Free evolution step size (unit: fs)
# }

# Molecule parameter dictionary in the document
# molecule_dictionary = {
#     'molecule_name': 'N2',
#     'T': 25,  # Temperature (unit: K)
#     'Jmax': 30,  # Maximum J value
#     'B': 1.9895,  # Rotational constant (unit: cm^-1)
#     'centrifugal': 0.0,  # Centrifugal distortion constant
#     'alphaper': 9.8,  # Perpendicular polarizability (unit: Bohr^3)
#     'alphapar': 15,  # Parallel polarizability (unit: Bohr^3)
#     'reneven': 0.6667,  # Weight for even J states
#     'renodd': 0.3333,  # Weight for odd J states
#     'tau': 80,  # First pulse width (unit: fs)
#     'I': 20,  # Peak intensity of the first pulse (unit: TW/cm^2)
#     'tau2': 50,  # Second pulse width (unit: fs)
#     'I02': 3,  # Peak intensity of the second pulse (unit: TW/cm^2)
#     'Tfree1': 9000,  # Free evolution time after the first pulse (unit: fs)
#     'Tfree2': 10,  # Free evolution time after the second pulse (unit: fs)
#     'tau3': 10,  # Third pulse width (unit: fs)
#     'Tfree3': 10,  # Free evolution time after the third pulse (unit: fs)
#     'I03': 0.0,  # Peak intensity of the third pulse (unit: W/cm^2)
#     'NOPlaser': 25,  # Number of time points in the laser field
#     'stepfree': 20  # Free evolution step size (unit: fs)
# }

# # Molecule parameter dictionary in the file (N2)
# molecule_dictionary = {
#     'molecule_name': 'n2_test',
#     'T': 298,  # Temperature (unit: K)
#     'Jmax': 65,  # Maximum J value
#     'B': 1.9895,  # Rotational constant (unit: cm^-1)
#     'centrifugal': 0.0,  # Centrifugal distortion constant
#     'alphaper': 9.8,  # Perpendicular polarizability (unit: Bohr^3)
#     'alphapar': 15,  # Parallel polarizability (unit: Bohr^3)
#     'reneven': 0.6667,  # Weight for even J states
#     'renodd': 0.3333,  # Weight for odd J states
#     'tau': 50,  # First pulse width (unit: fs)
#     'I': 3,  # Peak intensity of the first pulse (unit: TW/cm^2)
#     'tau2': 50,  # Second pulse width (unit: fs)
#     'I02': 3,  # Peak intensity of the second pulse (unit: TW/cm^2)
#     'Tfree1': 9000,  # Free evolution time after the first pulse (unit: fs)
#     'Tfree2': 10,  # Free evolution time after the second pulse (unit: fs)
#     'tau3': 10,  # Third pulse width (unit: fs)
#     'Tfree3': 10,  # Free evolution time after the third pulse (unit: fs)
#     'I03': 0.0,  # Peak intensity of the third pulse (unit: W/cm^2)
#     'NOPlaser': 25,  # Number of time points in the laser field
#     'stepfree': 20  # Free evolution step size (unit: fs)
# }

# # Molecule parameter dictionary in the file (CO2)
molecule_dictionary = {
    'molecule_name': 'co2_test',
    'T': 295,  # Temperature (unit: K)
    'Jmax': 65,  # Maximum J value
    'B': 0.39021,  # Rotational constant (unit: cm^-1)
    'centrifugal': 0.0,  # Centrifugal distortion constant
    'alphaper': 13.1,  # Perpendicular polarizability (unit: Bohr^3)
    'alphapar': 27.25,  # Parallel polarizability (unit: Bohr^3)
    'reneven': 1,  # Weight for even J states
    'renodd': 0,  # Weight for odd J states
    'tau': 50,  # First pulse width (unit: fs)
    'I': 3,  # Peak intensity of the first pulse (unit: TW/cm^2)
    'tau2': 50,  # Second pulse width (unit: fs)
    'I02': 3,  # Peak intensity of the second pulse (unit: TW/cm^2)
    'Tfree1': 80000,  # Free evolution time after the first pulse (unit: fs)
    'Tfree2': 10,  # Free evolution time after the second pulse (unit: fs)
    'tau3': 10,  # Third pulse width (unit: fs)
    'Tfree3': 10,  # Free evolution time after the third pulse (unit: fs)
    'I03': 1.0,  # Peak intensity of the third pulse (unit: W/cm^2)
    'NOPlaser': 25,  # Number of time points in the laser field
    'stepfree': 20  # Free evolution step size (unit: fs)
}

# Run the main function
cos2theta_original, D2, D4, D6, D8, time1, Eout, symm, mol_name = alignlinearD(molecule_dictionary)
cos2theta, D2, D4, D6, D8, time2, Eout, symm, mol_name = alignlinearD_DSIC(molecule_dictionary)

# Plot <cos^2 theta> vs time
plt.figure(figsize=(10, 6))
plt.plot(time1, cos2theta_original, label='Without DSIC', color='blue')  # Before correction
plt.plot(time2, cos2theta, label='With DSIC', color='red')  # After correction
plt.xlabel('Time (ps)')
plt.ylabel('<cos^2 theta>')
plt.title('Molecular Alignment with Dipole Self-Interaction Correction')
plt.legend()
plt.grid(True)

# Save the image as a PNG file
plt.savefig(f'../result/{mol_name}_cos2theta_vs_time_DSIC.png', dpi=300, bbox_inches='tight')