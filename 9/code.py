import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

D = np.array([228, 238, 248, 258, 268, 278, 280, 288, 260, 262, 264, 266, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 280, 283])
I = np.array([6.6, 5.9, 2.5, 1.7, 0.3, 0.2, 0.6, 6.2, 1.2, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1, 0.1, 0.0, 0.0, 0.1, 0.1, 0.2, 0.3, 0.5, 1.6])

D = D - 218

# Uncertainties
D_err = math.sqrt(2)
I_err = 0.05
sigma = np.full(D.shape, I_err, dtype=float)





font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)


# Define the model function with an amplitude scaling parameter A.
def mxc(x, n1, n2, A):
    # Convert angle to radians
    x_rad = x * math.pi / 180
    term = n1 / n2 * np.sin(x_rad)
    sqrt_term = np.sqrt(np.maximum(0, 1 - term**2))
    base = ((n1 * sqrt_term - n2 * np.cos(x_rad)) / (n1 * sqrt_term + n2 * np.cos(x_rad)))**2
    return A * base


# Using the gotten 0 values: 

print(f"REflective index is  {np.tan((273.5 - 218)*math.pi/180)} +- {math.sqrt(3)*0.5*(np.cos((273.5 - 218)*math.pi/180))**-2*math.pi/180} times that of air")








# Initial guess: n1, n2, and amplitude A
initial_guess = [1, 1.5, 1]

# Set bounds: n1 and n2 in [1,2], amplitude A in [0, infinity)
bounds = ([1, 1, 0], [2, 2, np.inf])

# Perform the curve fitting using curve_fit
popt, pcov = curve_fit(mxc, D, I, p0=initial_guess, sigma=sigma, bounds=bounds, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))

# Print the optimized parameters and their uncertainties
print("Optimized Parameters:")
print("n1 =", popt[0], "+/-", perr[0])
print("n2 =", popt[1], "+/-", perr[1])
print("Amplitude A =", popt[2], "+/-", perr[2])

# Plot the data and the best-fit curve
plt.errorbar(D, I, xerr=D_err, yerr=I_err, fmt='o', color="black", capsize=1, label="Data Points")
x_fit = np.linspace(0, max(D) + 2, 1000)
y_fit = mxc(x_fit, *popt)
plt.plot(x_fit, y_fit, label='Best Fit', color='red')
plt.xlabel('Angle (Degree)')
plt.ylabel('Current (mA)')
plt.legend()
plt.show()

























