Since `alignlinearD_CT.py` and `alignlinearmultipulses.py` are versions without the correction terms, all functions and parameters are essentially equivalent to those in `alignlinearD_CT_DSIC.py` and `alignlinearmultipulses_DSIC.py`. Therefore, the following explanation of functions will only focus on the latter two.

---

## **File `alignlinearD_CT_DSIC.py` Function Explanations**

### **Function `alignlinearD_DSIC(molecule_dictionary)`**

#### **Function**
Calculates the alignment behavior of molecules under the influence of laser pulses, incorporating dipole self-interaction corrections (DSIC).

#### **The core function**
1. Extracts parameters such as rotational constants, temperature, polarizability, etc., from the molecular configuration dictionary.
2. Calculates the electric field intensity of the laser pulse and the induced potential energy.
3. Calls the `alignlinearmultipulses` function to compute the alignment behavior of the molecules.
4. Plots the time evolution of the molecular alignment expectation value $ \langle \cos^2 \theta \rangle $ .

---

## **Function Explanations for `alignlinearmultipulses_DSIC.py`**

### **Function \( V(l_1, l_2, m) \)**

#### **Purpose**
Calculates the matrix element \( \langle Y_{l_1m} | \cos^2 \theta | Y_{l_2m} \rangle \), which represents the \( \cos^2 \theta \) matrix element between spherical harmonic functions.

#### **Core Role**
Used to construct the \( \cos^2 \theta \) matrix, describing the coupling between molecular rotational states.

---

### **Function \( V_6(l_1, l_2, m) \)**

#### **Purpose**
Calculates the matrix element \( \langle Y_{l_1m} | \cos^6 \theta | Y_{l_2m} \rangle \).

#### **Core Role**
Used to construct the \( \cos^6 \theta \) matrix, describing higher-order rotational state couplings.

---

### **Function \( V_8(l_1, l_2, m) \)**

#### **Purpose**
Calculates the matrix element \( \langle Y_{l_1m} | \cos^8 \theta | Y_{l_2m} \rangle \).

#### **Core Role**
Used to construct the \( \cos^8 \theta \) matrix, describing higher-order rotational state couplings.

---

### **Function \( \text{fieldode}(t_{\text{field}}, \text{coeff}, t, n, V_{00}, V_{02}, V_{04}, E_j, V_{\cos^2\theta}, V_{\cos^4\theta}) \)**

#### **Purpose**
Defines the field function for solving differential equations, considering the rotational Hamiltonian \( H_{\text{rot}} \), induced Hamiltonian \( H_{\text{ind}} \), and dipole self-interaction Hamiltonian \( H_{\text{self}} \).

#### **Core Role**
1. Constructs the Hamiltonian incorporating rotational energy, induced potential energy, and dipole self-interaction potential energy.
2. Solves the time-dependent Schrödinger equation to compute the time evolution of molecular state coefficients.

---

### **Function \( \text{alignlinearmultipulses}(J_{\text{max}}, J_0, m_0, \text{rotcon}, \text{centrifugal}, V_{00}, V_{02}, V_{04}, t, t_{\text{free}}, \text{NOPlaser}, \text{NOPfree}, \text{fullsize}, \text{pulsed}) \)**

#### **Purpose**
Calculates the alignment behavior of molecules under multiple laser pulses, considering dipole self-interaction corrections.

#### **Core Role**
1. Constructs coupling matrices between rotational states (e.g., \( \cos^2 \theta \), \( \cos^4 \theta \), etc.).
2. Calls the `fieldode` function to solve for the time evolution of molecular state coefficients.
3. Computes the alignment expectation values, such as \( \langle \cos^2 \theta \rangle \) and \( \langle \cos^4 \theta \rangle \).

---

## **Summary**

- **`alignlinearD_DSIC`**: The main function responsible for parameter initialization, laser field setup, induced potential energy calculation, and invoking other functions to compute alignment behavior.
- **`V`, `V6`, `V8`, `Vprime`**: Calculate the matrix elements of \( \cos^n \theta \) for different orders, used to describe the coupling between molecular rotational states.
- **`fieldode`**: Defines the differential equations and solves for the time evolution of molecular state coefficients.
- **`alignlinearmultipulses`**: Computes molecular alignment behavior under multiple laser pulses, incorporating dipole self-interaction corrections.

Using these functions, the program simulates the alignment behavior of molecules under laser pulses and considers dipole self-interaction corrections, providing a more accurate description of experimental phenomena.

## **Appendix: What Was Changed**

The following details explain how the formulas for \( V_{00} \), \( V_{02} \), and \( V_{04} \) were derived from the theoretical framework in the documentation.

---

### **1. Expansion of the Induced Potential Energy \( H_{\text{ind}} \)**

In the documentation, the expression for the induced potential energy \( H_{\text{ind}} \) is:

$$
H_{\text{ind}}(t) = \left( \Delta \alpha \cos^2 \theta + \alpha_{\perp} \right) \cdot \frac{\epsilon^2(t)}{4} \hat{Z}
$$

Where:
- \( \Delta \alpha = \alpha_{\parallel} - \alpha_{\perp} \)
- \( \alpha_{\perp} \) is the perpendicular polarizability.
- \( \epsilon(t) \) represents the time-dependent external electric field intensity.

#### **Expanding \( H_{\text{ind}} \)**

The induced potential \( H_{\text{ind}} \) can be expanded into two terms:

$$
H_{\text{ind}}(t) = \left( \Delta \alpha \cos^2 \theta + \alpha_{\perp} \right) \cdot \frac{\epsilon^2(t)}{4}
$$

This can be decomposed as:

$$
H_{\text{ind}}(t) = V_{00} + V_{02} \cos^2 \theta
$$

Where:
- \( V_{00} = \alpha_{\perp} \cdot \frac{\epsilon^2(t)}{4} \)
- \( V_{02} = \Delta \alpha \cdot \frac{\epsilon^2(t)}{4} \)

Thus:

$$
V_{00} = -\frac{1}{4} \alpha_{\perp} \epsilon^2(t)
$$
$$
V_{02} = -\frac{1}{4} (\alpha_{\parallel} - \alpha_{\perp}) \epsilon^2(t)
$$

---

### **2. Expansion of the Dipole Self-Interaction Correction \( H_{\text{self}} \)**

In the documentation, the expression for the dipole self-interaction correction \( H_{\text{self}} \) is:

$$
H_{\text{self}} = \lambda^2 \mu_Z^2 = \lambda^2 \left( \Delta \alpha \cos^2 \theta + \alpha_{\perp} \right)^2 \cdot \frac{\epsilon^2(t)}{4}
$$

Where:
- \( \lambda \) is the coupling constant for the dipole self-interaction correction.
- \( \mu_Z = \Delta \alpha \cos^2 \theta + \alpha_{\perp} \) is the induced dipole moment.

#### **Expanding \( H_{\text{self}} \)**

Expanding \( H_{\text{self}} \) into a polynomial:

$$
H_{\text{self}} = \lambda^2 \left( \Delta \alpha^2 \cos^4 \theta + 2 \alpha_{\perp} \Delta \alpha \cos^2 \theta + \alpha_{\perp}^2 \right) \cdot \frac{\epsilon^2(t)}{4}
$$

This can be decomposed as:

$$
H_{\text{self}} = V_{04} \cos^4 \theta + V_{06} \cos^2 \theta + V_{08}
$$

Where:
- \( V_{04} = \lambda^2 \Delta \alpha^2 \cdot \frac{\epsilon^2(t)}{4} \)
- \( V_{06} = \lambda^2 2 \alpha_{\perp} \Delta \alpha \cdot \frac{\epsilon^2(t)}{4} \)
- \( V_{08} = \lambda^2 \alpha_{\perp}^2 \cdot \frac{\epsilon^2(t)}{4} \)

Thus:

$$
V_{04} = -\frac{1}{4} \lambda^2 (\alpha_{\parallel} - \alpha_{\perp})^2 \epsilon^2(t)
$$
$$
V_{06} = -\frac{1}{4} \lambda^2 2 \alpha_{\perp} (\alpha_{\parallel} - \alpha_{\perp}) \epsilon^2(t)
$$
$$
V_{08} = -\frac{1}{4} \lambda^2 \alpha_{\perp}^2 \epsilon^2(t)
$$

---

### **3. Summary of Formulas**

| Induced Potential Term | Formula |
|-------------------------|---------|
| \( V_{00} \)           | \( -\frac{1}{4} \alpha_{\perp} \epsilon^2(t) \) |
| \( V_{02} \)           | \( -\frac{1}{4} (\alpha_{\parallel} - \alpha_{\perp}) \epsilon^2(t) \) |
| \( V_{04} \)           | \( -\frac{1}{4} \lambda^2 (\alpha_{\parallel} - \alpha_{\perp})^2 \epsilon^2(t) \) |
| \( V_{06} \)           | \( -\frac{1}{4} \lambda^2 2 \alpha_{\perp} (\alpha_{\parallel} - \alpha_{\perp}) \epsilon^2(t) \) |
| \( V_{08} \)           | \( -\frac{1}{4} \lambda^2 \alpha_{\perp}^2 \epsilon^2(t) \) |

---

### **4. Summary**

- \( V_{00} \) and \( V_{02} \) are interaction terms between the induced dipole moment and the external electric field.
- \( V_{04} \), \( V_{06} \), and \( V_{08} \) are dipole self-interaction correction terms.
- These formulas are derived by expanding the expressions for the induced potential energy and dipole self-interaction corrections.
- In the code, \( V_{06} \) and \( V_{08} \) are of the same order as \( V_{02} \) and \( V_{00} \), respectively, so like terms are combined.
