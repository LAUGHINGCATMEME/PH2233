import math
import numpy as np
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
   

