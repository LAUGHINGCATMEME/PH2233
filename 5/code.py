import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import RealData, Model, ODR

l = np.array([16, 18, 20, 24, 30, 35, 40])
i = np.array([12.2, 10.0, 7.6, 4.9, 3.3, 2.5, 2.2])
l = l * 10**-2
i = i * 10**-6
ln_l = np.log(l)
ln_i = np.log(i)
l_least = 0.001
i_least = 10**-7
ln_l_err = (l_least /2)/l
ln_i_err = (i_least /2)/i




def mxc(params, x):
        m, c  = params
        return m * x + c

odr_instance = ODR(RealData(ln_l, ln_i, sx=ln_l_err, sy=ln_i_err), Model(mxc), beta0=[-1, 1])
out = odr_instance.run()
print("Optimized Parameters:", out.beta)
print("Parameter Errors:", out.sd_beta)
plt.errorbar(ln_l, ln_i, marker='o', linestyle='', xerr=ln_l_err, yerr=ln_i_err, color="black", ecolor='black', capsize=1, label="Data Points")
x_fit = np.linspace(min(ln_l), max(ln_l), 1000)
y_fit = mxc(out.beta, x_fit)
plt.plot(x_fit, y_fit, label='Linear ODR Fir')
plt.xlabel('ln(Distance (m))')
plt.ylabel('ln(Current (A))')
plt.legend()






# Provided data (wavelengths in nm and current in A)
l = np.array([635, 570, 540, 500, 460])
i = np.array([0.32, 0.48, 0.65, 0.71, 0.99])

# Convert wavelength to meters (SI units)
l = l * 1e-9

# Measurement error estimates
l_least = 1e-8
i_least = 0.01
l_err = l_least / 2
i_err = i_least / 2

c = 299792458  # Speed of light in m/s

# Compute frequency and its error: v = c/l
v = c / l
v_err = c * l_err / l**2


# Optional: Scale frequency to improve numerical stability (e.g., divide by 1e14)
v_scaled = v / 1e14
v_err_scaled = v_err / 1e14

# Define the linear model with two parameters (slope and intercept)
def linear_model(params, x):
    m, b = params
    return m * x + b

# Set up ODR with the scaled frequency data
data = RealData(v_scaled, i, sx=v_err_scaled, sy=i_err)
model = Model(linear_model)
odr_instance = ODR(data, model, beta0=[-1, 1])
out = odr_instance.run()

print("Optimized Parameters:", out.beta)
print("Parameter Errors:", out.sd_beta)

# Plot the data and the fitted line
plt.errorbar(v_scaled, i, xerr=v_err_scaled, yerr=i_err, fmt='o', 
             color='black', ecolor='black', capsize=1, label="Data Points")
x_fit = np.linspace(0, max(v_scaled), 1000)
y_fit = linear_model(out.beta, x_fit)
plt.plot(x_fit, y_fit, label='Linear ODR Fit')
plt.xlabel('Frequency (v / 1e14 Hz)')
plt.ylabel('Current (A)')
plt.legend()
plt.grid()
plt.show()



