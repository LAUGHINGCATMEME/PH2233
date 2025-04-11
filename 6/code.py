import numpy as np 
import pandas as pd
import math
import csv

# VARIABLES
d = 5#mm
p = 
pa = 
T = 
P = 
n = 
c = 
g = 

#############DATA FRAME ######################
dynamic = pd.read_csv('dynamic.csv')
balancing = pd.read_csv('balancing.csv')

# CONVERT TO SI
d = d * 10**-3


# TABLE FILL
d_tf = 
d_tr = 
d_vf = 
d_V = 

b_tf = 
b_vf = 
b_


# OTHER VARIABLES CALC
C = 4 * math.pi * d * g * (p - pa) / 3 
D = 9 * n / (2 * g * (p - pa)) 
Z = c/(2 * P)



print(dynamic)


