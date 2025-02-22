''' Directional Microwave Catheter Antenna '''
''' Arav Sharma '''
''' Sarjit Bharj '''
''' Princeton Microwave Technology '''

import numpy as np
import matplotlib.pyplot as plt

def elliptic_integral_ratio(k):
    k_prime = np.sqrt(1 - k**2)
    ratio = np.pi / np.log(2 * (1 + np.sqrt(k_prime)) / (1 - np.sqrt(k_prime)))
    return ratio

def calculate_characteristic_impedance(a, b, c, S, W, eps1, eps2, eps_in, eps_out):
    eps0 = 8.854e-12  # permittivity of free space

    # calculate k_a
    k_a = (S / (S + 2 * W)) * np.sqrt((1 - (S + 2 * W)**2/ (4 * b**2 * np.pi**2)) / (1 - S**2 / (4 * b**2 * np.pi**2)))

    # Calculate elliptic integral ratio for k_a
    K_ratio_a = elliptic_integral_ratio(k_a)

    # Calculate C^a
    C_a = 4 * eps0 * K_ratio_a

    # Calculate effective permittivity contributions
    A1 = np.pi / (4 * b * np.log(b / a))
    k_1 = (np.sinh(A1 * S) / np.sinh(A1 * (S + 2 * W))) * np.sqrt(
        (1 - (np.sinh(A1 * (S + 2 * W))**2 / np.sinh(2 * A1 * b * np.pi)**2)) /
        (1 - (np.sinh(A1 * S)**2 / np.sinh(2 * A1 * b * np.pi)**2))
    )
    K_ratio_1 = elliptic_integral_ratio(k_1)
    C_s1_a = 2 * eps0 * K_ratio_1

    A2 = np.pi / (4 * b * np.log(b / c))
    k_2 = (np.sinh(A2 * S) / np.sinh(A2 * (S + 2 * W))) * np.sqrt(
        (1 - (np.sinh(A2 * (S + 2 * W))**2 / np.sinh(2 * A2 * b * np.pi)**2)) /
        (1 - (np.sinh(A2 * S)**2 / np.sinh(2 * A2 * b * np.pi)**2))
    )
    K_ratio_2 = elliptic_integral_ratio(k_2)
    C_s2_a = 2 * eps0 * K_ratio_2

    # Calculate q_in and q_out
    q_in = (C_a / 2 - C_s1_a) / C_a
    q_out = (C_a / 2 - C_s2_a) / C_a

    # Calculate effective permittivity
    eps_eff = q_in * eps_in + q_out * eps_out

    # Calculate characteristic impedance
    Z_c = (120 * np.pi * eps0) / (C_a * np.sqrt(eps_eff))

    return Z_c

# Define parameters
a = 10
b = 35
c = 42
S = 10
eps1 = 2.2
eps2 = 2.2
eps_in = 1
eps_out = 47 + 16j

W = np.linspace(10, 70, 50)

Z_c = []

for w in W:
    Z_c += [calculate_characteristic_impedance(a, b, c, S, w, eps1, eps2, eps_in, eps_out)]

print(Z_c)

plt.figure()

plt.title('Zc vs. W')

plt.ylabel('Zc (Ohms)')

plt.xlabel('W (mils)')

plt.plot(W, np.real(Z_c), 'r--')
plt.plot(W, np.imag(Z_c), 'b')

plt.grid()

plt.savefig('Z_cVsW.png')