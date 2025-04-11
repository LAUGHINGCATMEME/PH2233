#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt

# ================================================================
# Constants
# ================================================================
# Diffraction grating: 15000 lines per inch => spacing d = (1 inch)/(15000)
inch_in_m = 0.0254         # 1 inch in m
lines_per_inch = 15000
d = inch_in_m / lines_per_inch  # Grating spacing in m
m_order = 1                 # Use first order diffraction

# Accepted Rydberg constant
R_accepted = 1.097e7        # in m^-1

# Other constants for additional calculations:
h = 6.6260755e-34           # Planck constant in J.s
c = 3.0e8                   # Speed of light in m/s
alpha = 1/137.0             # Fine structure constant ~1/137
eV_in_J = 1.60e-19          # 1 eV in Joules

# ================================================================
# Helper Functions
# ================================================================
def deg_min_to_deg(deg, minute):
    """Convert degrees and minutes to a floating point degree value."""
    return deg + minute/60.0

def angle_in_radians_from_deg_min(deg, minute):
    """Convert (deg, minute) tuple into radians."""
    return math.radians(deg_min_to_deg(deg, minute))

def effective_angle(measured, zero):
    """
    Given measured and zero values as (deg, minute),
    compute the difference (in degrees).
    """
    measured_deg = deg_min_to_deg(*measured)
    zero_deg = deg_min_to_deg(*zero)
    return measured_deg - zero_deg

# ================================================================
# Mercury Lamp Data: (deg, min) readings
# The zero (reference) for mercury lamp is given as (186,22)
# ================================================================
zero_mercury = (186, 22)
mercury_data = [
    {"wavelength_nm": 405, "reading": (186, 22)},
    {"wavelength_nm": 436, "reading": (198, 8)},
    {"wavelength_nm": 546, "reading": (199, 2)},
    {"wavelength_nm": 577, "reading": (201, 47)},
    {"wavelength_nm": 579, "reading": (202, 49)}
]

print("=== Mercury Lamp Data ===")
mercury_results = []
for entry in mercury_data:
    # Compute effective angle in degrees then in radians
    eff_deg = effective_angle(entry["reading"], zero_mercury)
    eff_rad = math.radians(eff_deg)
    # Use diffraction grating equation: lambda = d*sin(theta) for m = 1
    lambda_exp = d * math.sin(eff_rad)
    # Convert measured wavelength from nm for comparison (if desired)
    lambda_measured_nm = entry["wavelength_nm"]
    mercury_results.append((eff_deg, lambda_exp, lambda_measured_nm))
    print(f"Wavelength {lambda_measured_nm:3d} nm: effective angle = {eff_deg:6.2f}°; "
          f"Computed lambda = {lambda_exp*1e9:6.2f} nm")
print()

# ================================================================
# Hydrogen Lamp Data:
# Initial (zero) reading is given as 264° (assumed 264°00')
# For each color we have a clockwise reading and an anticlockwise reading.
# Effective angle is computed as:
#    For clockwise: (clockwise - initial) ; for anticlockwise: (initial - anticlockwise)
# Then take the average.
# ================================================================
zero_hydrogen = (264, 0)
hydrogen_data = [
    {"color": "Red",    "cw": (282, 24), "acw": (244, 23), "nf": 3},
    {"color": "Yellow", "cw": (280, 56), "acw": (246, 27), "nf": 4},
    {"color": "Green",  "cw": (277, 25), "acw": (249, 35), "nf": 5},
    {"color": "Violet", "cw": (276, 49), "acw": (250, 15), "nf": 6}
]

print("=== Hydrogen Lamp Data ===")
hydro_angles = []   # to store effective angles for regression
hydro_lambdas = []  # to store computed lambda values from hydrogen data

for entry in hydrogen_data:
    # Calculate effective angle in degrees (clockwise and anticlockwise)
    # Clockwise: measured - initial; Anticlockwise: initial - measured.
    eff_cw_deg = deg_min_to_deg(*entry["cw"]) - deg_min_to_deg(*zero_hydrogen)
    eff_acw_deg = deg_min_to_deg(*zero_hydrogen) - deg_min_to_deg(*entry["acw"])
    eff_deg_avg = (eff_cw_deg + eff_acw_deg) / 2.0
    eff_rad = math.radians(eff_deg_avg)
    # Compute lambda using diffraction equation
    lambda_exp = d * math.sin(eff_rad)
    hydro_angles.append(eff_deg_avg)
    hydro_lambdas.append(lambda_exp)
    print(f"{entry['color']:<7s} (nf = {entry['nf']}): cw = {eff_cw_deg:5.2f}°, "
          f"acw = {eff_acw_deg:5.2f}°, avg = {eff_deg_avg:5.2f}°; "
          f"Computed lambda = {lambda_exp*1e9:6.2f} nm")
print()

# ================================================================
# Linear Regression to Determine the Rydberg Constant
# For the Balmer series with n_i = 2, we have:
#    1/lambda = R (1/2^2 - 1/n_f^2)
# We form x = (1/4 - 1/n_f^2) and y = 1/lambda for each data point.
# Then, R_exp = slope obtained from regression forced through the origin.
# ================================================================
n_initial = 2
x_vals = []
y_vals = []
for entry, lam in zip(hydrogen_data, hydro_lambdas):
    n_final = entry["nf"]
    x = (1/ (n_initial**2)) - (1/(n_final**2))  # note: 1/4 - 1/n_f^2
    y = 1/lam   # lambda in m so y in m^-1
    x_vals.append(x)
    y_vals.append(y)

x_vals = np.array(x_vals)
y_vals = np.array(y_vals)

# Perform a least-squares fit forced through the origin.
# The slope is given by: slope = sum(x*y)/sum(x^2)
slope = np.sum(x_vals*y_vals) / np.sum(x_vals**2)
R_measured = slope

print("=== Determination of the Rydberg Constant ===")
print("Data (for regression):")
print(" nf |    x = (1/4 - 1/nf^2)    | 1/lambda (m^-1)")
for entry, x_val, y_val in zip(hydrogen_data, x_vals, y_vals):
    print(f"{entry['nf']:2d} | {x_val:20.6f} | {y_val:12.2e}")

print(f"\nExperimental Rydberg Constant: R_exp = {R_measured:.3e} m^-1")
percent_diff_R = abs(R_accepted - R_measured)/R_accepted * 100
print(f"Percent difference from accepted value: {percent_diff_R:.2f} %\n")

# ================================================================
# Compare Experimental Wavelengths vs Predicted from Measured R
# For each Balmer line:
#    lambda_pred = 1 / [R_measured * (1/4 - 1/nf^2)]
# ================================================================
print("=== Wavelength Comparison (Hydrogen) ===")
print("Color   nf   Computed lambda (nm)   lambda_pred (nm)   % difference")
for entry, lam in zip(hydrogen_data, hydro_lambdas):
    n_final = entry["nf"]
    x = (1/ (n_initial**2)) - (1/(n_final**2))
    lambda_pred = 1/(R_measured * x)
    pct_diff = abs(lambda_pred - lam)/lambda_pred * 100
    print(f"{entry['color']:<7s} {n_final:2d}   {lam*1e9:18.2f}    {lambda_pred*1e9:14.2f}    {pct_diff:8.2f} %")
print()

# ================================================================
# Additional Calculations
# ================================================================
# 1. Bohr Radius:
# Using: r1 = 4*pi*epsilon0*hbar^2/(m_e*e^2)
# For typical values: r1 ≈ 5.29e-11 m
r1 = 5.29e-11  # (m)
print(f"Bohr Radius: r1 ≈ {r1:.3e} m")

# 2. Electron speed at the Bohr radius: v1 = alpha * c
v1 = alpha * c
print(f"Electron speed at r1: v1 ≈ {v1:.3e} m/s")

# 3. Photon Energies: for red photon (600 nm) and blue photon (400 nm)
lambda_red = 600e-9  # m
lambda_blue = 400e-9  # m
E_red = h * c / lambda_red      # in Joules
E_blue = h * c / lambda_blue    # in Joules
E_red_eV = E_red / eV_in_J
E_blue_eV = E_blue / eV_in_J

print("\nPhoton Energies:")
print(f"Red photon (600 nm):  {E_red:.3e} J   ~ {E_red_eV:.2f} eV")
print(f"Blue photon (400 nm): {E_blue:.3e} J   ~ {E_blue_eV:.2f} eV")

# ================================================================
# OPTIONAL: Plot the regression for the Balmer series
# ================================================================
plt.figure(figsize=(6,4))
plt.plot(x_vals, y_vals, 'bo', label='Data points')
# Line through the origin with slope = R_measured:
x_fit = np.linspace(0, np.max(x_vals)*1.1, 100)
y_fit = R_measured * x_fit
plt.plot(x_fit, y_fit, 'r-', label=f'Fit: slope = {R_measured:.2e} m^-1')
plt.xlabel(r'$1/4 - 1/n_f^2$')
plt.ylabel(r'$1/\lambda$ (m$^{-1}$)')
plt.title("Linear Regression for Rydberg Constant Determination")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

