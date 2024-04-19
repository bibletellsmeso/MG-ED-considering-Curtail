# library import
import cvxpy as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%--------------------------Hawaai Problem with exact parameters----------#
# Generate data
P_PV = np.loadtxt("PV_for_scheduling.txt")               # 6/7 Solar data Input
P_L = np.loadtxt("Load_for_scheduling.txt")             # 6/7 Load data Input

T = len(P_L)                                            # Total time slots as 1 slot is 15 min so in total it will be 24*4
t = 0.25                                                # Delta t is 1/4 hr

# DG data
a = 32000
b = 210
c = 0.097
# P_DGmin = 225                                           # Diesel is always on, minimal power limit 30% of maximum
P_DGmax = 750
i = 0
piece = 10
P_DGpmax = P_DGmax/piece

def F_cost(p):
    return (c*(p**2) + b*p + a)

mincost = F_cost(0)

Sl_diesel = []                                    # slope of each piece
for i in range(piece):
    Sl_diesel.append((F_cost(P_DGpmax*(i+1)) - F_cost(P_DGpmax*i)) / P_DGpmax)
#%%--------------------------Battery parameters----------#
ESS_max = 250
ESS_cap = 750
C_eff = 93/100
D_eff = 100/93

#%% Variables intialization
# Unit
P_DG = cp.Variable(T)                                  # DG variable
u = cp.Variable(T, integer=True)                       # ON/OFF
r = []                                                 # "output k-th segment at time step i" r[]][]:(10,96) /
for i in range(10):
    r.append(cp.Variable(T))

# P_Curt = cp.Variable(T)                                # PV Curtailment

      

# ESS 
ESS_char = cp.Variable(T)
ESS_dischar = cp.Variable(T)
ESS = cp.Variable(T)

Uc = cp.Variable(T, integer=True)
Ud = cp.Variable(T, integer=True)

SOC = cp.Variable(T)

#%% constraints
constraints = []
for i in range(T):
    temp = []
    # Load balance
    temp.append(P_DG[i] + ESS[i] + P_PV[i] == P_L[i]) # - P_Curt[i]
    
    # ESS
    temp.append(Uc[i] >= 0) 
    temp.append(Uc[i] <= 1)
    temp.append(Ud[i] >= 0)
    temp.append(Ud[i] <= 1)
    temp.append(Uc[i] + Ud[i] <= 1)
            
    temp.append(ESS_char[i] >= -1.0 * (ESS_max * Uc[i]))
    temp.append(ESS_char[i] <= 0)
    temp.append(ESS_dischar[i] <= ESS_max * Ud[i])
    temp.append(ESS_dischar[i] >= 0)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    
    temp.append(ESS[i] == ESS_char[i] + ESS_dischar[i])
    
    # SOC
    temp.append(SOC[i] >= 20)
    temp.append(SOC[i] <= 80)
    temp.append(SOC[i] == SOC[i-1] - ((ESS_char[i] * 100 * C_eff * t) / (ESS_cap)) - ((ESS_dischar[i] * 100 * D_eff * t) / (ESS_cap))
                if i >= 1 else SOC[i] == 50)
    temp.append(SOC[95] == 50)

        # segment
    for k in range(piece):
        temp.append(r[k][i] >= 0)
        temp.append(r[k][i] <= P_DGpmax * u[i])
    
    # Generation capacity
    temp.append(P_DG[i] >= 0)                     # 225, always ON
    temp.append(P_DG[i] <= 750 * u[i])            # Max
    temp.append(P_DG[i] == r[0][i] + r[1][i] + r[2][i] + r[3][i] + r[4][i] + r[5][i] + r[6][i] + r[7][i] + r[8][i] + r[9][i])
    
        # integer
    temp.append(u[i] >= 0)
    temp.append(u[i] <= 1)

        
    # PV curtail
    # temp.append(P_Curt[i] = c1[i])
        
    # sum every time
    constraints += temp
#%% objective function
cost = 0

obj = 0
for i in range(T):
    obj += (mincost * u[i]) * t
    for j in range(10):
        obj += (Sl_diesel[j] * r[j][i]) * t
#%% Optimal Schedule
prob = cp.Problem(cp.Minimize(obj), constraints)
prob.solve(solver=cp.GUROBI)
optimal_cost = np.round(prob.value, 3)
print('cost :', optimal_cost, 'KRW')
print(prob.status)
#%% Data for Plot
SOC_val = SOC.value
P_DG_val = P_DG.value
ESS_char_val = ESS_char.value
ESS_dischar_val = ESS_dischar.value
ESS_val = ESS.value
# print(p.value)
time = list(range(1, 97))

#%% Plot
Fig = plt.figure(figsize=[15,10])
ax1 = plt.subplot(2, 1, 2)
plt.xlabel('Time (h)')
plt.ylabel('OutPut (kW)')
plt.suptitle("ESS Scheduling", fontsize=25)
plt.plot(time, P_L, 'g-', label = 'Load')
plt.plot(time, P_DG_val, 'r-', label = 'Diesel')
plt.plot(time, ESS_val, 'b-', label = 'ESS')
plt.plot(time, P_PV, 'y-', label = 'PV')
plt.legend(loc=2)

plt.subplot(2, 1, 1, sharex=ax1)
plt.ylabel('SOC(%)')
plt.plot(SOC_val, 'b', label = 'SOC')
plt.xticks(visible=False)
plt.legend(loc=2)

plt.show()