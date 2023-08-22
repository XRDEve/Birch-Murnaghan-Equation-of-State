import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# Data preparation---------------------------------------------------------------------------------------------------------------
np.random.seed(9845)  # For reproducibility; controls random number generation.

# Requirements-------------------------------------------------------------------------------------------------------------------
P_exp, a_axis, a_err, b_axis, b_err, c_axis, c_err, V_exp, V_err = [], [], [], [], [], [], [], [], []

# Input Data from Experimental File----------------------------------------------------------------------------------------------
with open("PV.txt", "r") as XRData:   # input the name of the txt file. make sure it's in the correct directory
    linesXRData = XRData.readlines()
    for x in linesXRData:
        values = list(map(float, x.split('\t')))
        P_exp.append(values[0])
        a_axis.append(values[1])
        a_err.append(values[2])
        b_axis.append(values[3])
        b_err.append(values[4])
        c_axis.append(values[5])
        c_err.append(values[6])
        V_exp.append(values[7])
        V_err.append(values[8])

V_pred = np.arange(min(V_exp), max(V_exp) + 10)
V_pred = V_pred.reshape(-1, 1)

# SOLVE MODEL FOR THE 2ND-ORDER BMEoS
x0 = np.array([150, 1600])
xl = np.array([4, 0])
xu = np.array([1500, 1800])
options = {'algorithm': 'trust-region-reflective',
           'maxiter': 3000,
           'maxfev': 3000,
           'epsfcn': 1e-30,
           'xtol': 1e-30,
           'ftol': 1e-30,
           'disp': True}

# Create function for 2nd-order BMEoS
def fun2BMEoS(x):
    return (((((3/2)*x[0]) * (((x[1] / V_exp) ** (7/3)) - ((x[1] / V_exp) ** (5/3))) - P_exp) ** 2) / V_err)

# 'lm' for:  Levenberg–Marquardt algorithm. If bounds are required then use 'trf': trust-region-reflective  and  bounds=(xl, xu)
res = least_squares(objective, x,  bounds=(xl, xu), method='trf')
fitted_x = res.x

# Calculated Pressures for finer volume steps
P2BMEoS_data = ((3/2)*x[0]) * (((x[1] / V_exp) ** (7/3)) - ((x[1] / V_exp) ** (5/3)))
P2BMEoS = ((3/2)*x[0]) * (((x[1] / V_pred) ** (7/3)) - ((x[1] / V_pred) ** (5/3)))

SR2 = (P_exp - P2BMEoS_data) ** 2
SSR2 = np.sum(SR2)
R2_2 = np.linalg.lstsq(P2BMEoS_data, P_exp, rcond=None)[0]

# SOLVE MODEL FOR THE 3RD ORDER BMEoS
m0 = np.array([150, 4, 1500])
mlb = np.array([1, 1, 0])
mub = np.array([1500, 6, 1600])

