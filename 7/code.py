import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

T_a = 23.0

T_values = np.array([91.5, 90.5, 89.5, 88.5, 85.5, 84.5, 82.5, 81.5, 80.5, 79.5, 78.5, 77.5])
delta_t =   np.array([7,    8,    11,   11,   16,   17,   21,   21,   19,   23,   22,   22])

x = T_values - T_a

y = 1.0 / delta_t

def linear_zero_intercept(B, x):
    return B[0] * x

err_x = 0.25*math.sqrt(2)
err_y = 0.5 / delta_t**2

data = RealData(x, y, sx=err_x, sy=err_y)
model = Model(linear_zero_intercept)
beta0 = [0.02]  

odr_obj = ODR(data, model, beta0=beta0)
output = odr_obj.run()

k = output.beta[0]
k_err = output.sd_beta[0]

print(f"Best-fit cooling constant: k = {k:.6f} ± {k_err:.6f} s⁻¹")

T_target = 84.0
excess_target = T_target - T_a
cooling_rate = -k * excess_target
cooling_rate_err = excess_target * k_err  # error propagation

print(f"At T = {T_target}°C, the cooling rate dT/dt = {cooling_rate:.4f} ± {cooling_rate_err:.4f} °C/s")

# For visualization: Plot y (1/delta_t) vs x (T - T_a) and the best-fit line.
x_fit = np.linspace(np.min(x) - 1, np.max(x) + 1, 200)
y_fit = k * x_fit

plt.figure(figsize=(8,5))
plt.errorbar(x, y, yerr=err_y, fmt='o', label='Data: 1/Δt')
plt.plot(x_fit, y_fit, 'r-', label=f'Best fit: k = {k:.6f} s⁻¹')
plt.xlabel('Excess Temperature, (T - T_a) [°C]')
plt.ylabel('1/Δt [1/s]')
plt.title("Regression to Determine Cooling Constant k")
plt.legend()
plt.grid(True)
plt.show()







m = 0.905 
dm = 0.001
s = 0.380
ds = 0.001
Td = cooling_rate *-1 
dTd = cooling_rate_err

l = 331*10**-5
d = 9.90*10**-2
T2 = 95
T1 = T_target 


dd = 20**-4/2
dl = 10**-5/2
dT = 0.5/2
k = 4*m*s*Td/(math.pi*d**2*l*(T2-T1))
dk = k*math.sqrt((dm/m)**2+(ds/s)**2+(dTd/Td)**2+(2*dd/d)**2+(dl/l)**2+2*(dT/(T2-T1))**2)

print(f"{k} \pm {dk}")













exit()
# Original data list

data = [
    ("Aluminium", "24.0-23.0-24.7", 99.0, 59.9, 105),
    ("Brass",     "24.1-23.2-24.3", 99.0, 59.7, 85),
    ("Steel",     "22.1-24.8-20.5", 99.0, 60.0, 74),
    ("Aluminium", "24.3-23.7-24.3", 99.0, 59.8, 104),
    ("Brass",     "23.7-22.4-24.3", 99.0, 60.0, 85),
    ("Steel",     "24.6-25.3",      99.0, 59.9, 76),
    ("Brass",     "24.8-25.3",      99.0, 60.1, 86),
    ("Steel",     "23.3-23.5",      99.0, 59.7, 76),
    ("Copper",    "24.0",           99.0, 59.8, 75), 
    ("Copper",    "25.5-24.5-25.5", 99.0, 59.7, 74)
]

# Create a new data list by converting the second column of each entry into a list of floats.
data = [
    (material, float(np.average(list(map(float, range_str.split('-'))))), hundred, measurement, last)
    for material, range_str, hundred, measurement, last in data
]
delta_L_uncertainty = 1e-6 /2  # meters
L_uncertainty = 1e-3       # meters
temp_uncertainty = 1     /2  # Celsius

results = []

for row in data:
    material, t_initial, t_final, L_cm, delta_L_1e5 = row
    
    # Adjust the length 
    adjusted_L_cm = L_cm  
    adjusted_L = adjusted_L_cm * 1e-2  # Convert cm to meters
    
    # Calculate delta_T and handle division by zero if any
    delta_T = t_final - t_initial
    
    # Convert delta_L from 1e-5 m to meters
    delta_L = delta_L_1e5 * 1e-5
    
    # Compute expansion coefficient alpha
    alpha = delta_L / (adjusted_L * delta_T)
    
    # Calculate uncertainties
    sigma_delta_L = delta_L_uncertainty
    sigma_L = L_uncertainty
    sigma_delta_T = temp_uncertainty * math.sqrt(2)  # Combined uncertainty for delta_T
    
    rel_error_delta_L = sigma_delta_L / delta_L
    rel_error_L = sigma_L / adjusted_L
    rel_error_delta_T = sigma_delta_T / delta_T
    
    # Combined relative error for alpha
    rel_error_alpha = math.sqrt(rel_error_delta_L**2 + rel_error_L**2 + rel_error_delta_T**2)
    sigma_alpha = alpha * rel_error_alpha
    
    results.append((material, alpha, sigma_alpha))

for result in results:
    material, alpha, sigma_alpha = result
    print(f"Material: {material}")
    print(f"  α = {alpha:.4e} ± {sigma_alpha:.4e} 1/°C")
    print("-------------------")
   










