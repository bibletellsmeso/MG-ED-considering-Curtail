# library import
import cvxpy as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy

#%%--------------------------Hawaai Problem with exact parameters-------------#
# Generate data
P_PV = np.loadtxt("PV_for_scheduling.txt")                                     # 6/7 Solar data Input
P_L = np.loadtxt("Load_for_scheduling.txt")                                    # 6/7 Load data Input

T = len(P_L)                                                                   # Total time slots as 1 slot is 15 min so in total it will be 24*4
t = 15/60                                                                      # Delta t is 1/4 hr

# DG data
a = 32000
b = 210
c = 0.097
P_DGmin = 750*0.3                                                              # Diesel is always on, minimal power limit 30% of maximum
P_DGmax = 750
i = 0
piece = 10
P_DGpmax = P_DGmax/piece
TotalCurt = 600

def F_cost(p):
    return (c*(p**2) + b*p + a)                                                # 1/4 (KRW/hour) -> (KRW/min_15)

mincost = F_cost(0)

Sl_diesel = []                                                                 # slope of each piece
for i in range(piece):
    Sl_diesel.append((F_cost(P_DGpmax*(i+1)) - F_cost(P_DGpmax*i)) / P_DGpmax)
#%%--------------------------Battery parameters-------------------------------#
ESS_max = 500                                                                  # 250
ESS_CAP = 567                                                                  # 750
C_eff = 93/100
D_eff = 100/93
P_pump = 306/(24*4)                                                            # 306kW for 24 h
P_UPS = 4.5
#%%------------------------Variables intialization----------------------------#
# Unit
P_DG = cp.Variable(T)                                                          # DG variable
u = cp.Variable(T, integer=True)                                               # ON/OFF
p = []                                                                         # "output k-th segment at time step i" r[]][]:(10,96) /
for i in range(piece):
    p.append(cp.Variable(T))

P_Curt = cp.Variable(T)                                                        # PV Curtailment
P_PV_Output = cp.Variable(T)

# curt = cp.Variable(T)                                                        # 1 Piece

curt = []                                                                      # 10 Pieces
for i in range(piece):
    curt.append(cp.Variable(T))

# s1 = [1e-5]                                                                  # 1 Piece price
s10 = [1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5]            # 10 Pieces price

# ESS 
ESS_char = cp.Variable(T)
ESS_dischar = cp.Variable(T)
ESS = cp.Variable(T)

Uc = cp.Variable(T, integer=True)
Ud = cp.Variable(T, integer=True)

Upvc = cp.Variable(T, integer=True)

SOC = cp.Variable(T)

#%%----------------------------constraints------------------------------------#
constraints = []
for i in range(T):
    temp = []
    # Load balance
    temp.append(P_DG[i] + ESS[i] + P_PV[i] - P_Curt[i] == P_L[i])
    
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
    temp.append(SOC[i] == SOC[i-1] - ((ESS_char[i] * 100 * C_eff * t) / (ESS_CAP)) - ((ESS_dischar[i] * 100 * D_eff * t) / (ESS_CAP))
                if i >= 1 else SOC[i] == 50)
    temp.append(SOC[95] == 50)
    
    # Piecewise segment Power
    for j in range(piece):
        temp.append(p[j][i] >= 0)
        temp.append(p[j][i] <= P_DGpmax * u[i])
    
    # Generation capacity
    temp.append(u[i] >= 0)
    temp.append(u[i] <= 1)
    
    temp.append(P_DG[i] >= P_DGmin)                                            # always ON
    temp.append(P_DG[i] <= P_DGmax * u[i])
    temp.append(P_DG[i] == p[0][i] + p[1][i] + p[2][i] + p[3][i] + p[4][i] + p[5][i] + p[6][i] + p[7][i] + p[8][i] + p[9][i])

    # PV Curtailment
    temp.append(Upvc[i] >= 0)
    temp.append(Upvc[i] <= 1)

    temp.append(P_Curt[i] >= 0)
    temp.append(P_Curt[i] <= P_PV[i] * Upvc[i])

    # temp.append(curt[i] >= 0)                                                # 1 Piece
    # temp.append(curt[i] <= TotalCurt)
    
    for j in range(piece):                                                     # 10 Pieces
        temp.append(curt[j][i] >= 0)
        temp.append(curt[j][i] <= TotalCurt/piece)
    
    # # Piece Curtail 1 piece
    # temp.append(P_Curt[i] == curt[i])
    
    # Piece Curtail 10 pieces
    temp.append(P_Curt[i] == curt[0][i] + curt[1][i] + curt[2][i] + curt[3][i] + curt[4][i] + curt[5][i] + curt[6][i] + curt[7][i] + curt[8][i] + curt[9][i])
    temp.append(P_PV_Output[i] == P_PV[i] - P_Curt[i])
    
    # sum every time
    constraints += temp
#%%---------------------------------objective function------------------------#
cost = 0

obj = 0
for i in range(T):
    obj += t * mincost # * u[i]
    for j in range(piece):
        # obj += t * (Sl_diesel[j] * p[j][i] + s1 * curt[i])                   # 1 Piece
        obj += t * (Sl_diesel[j] * p[j][i] + s10[j] * curt[j][i])              # 10 Pieces
#%%----------------------------------Optimal Schedule-------------------------#
prob = cp.Problem(cp.Minimize(obj), constraints)
# prob.solve()
prob.solve(solver=cp.GUROBI)
optimal_cost = np.round(prob.value, 6)
print(P_DG.value)
print('cost :', optimal_cost, 'KRW')
print(prob.status)
#%% Data for Plot
SOC_val = SOC.value
P_DG_val = P_DG.value
ESS_char_val = ESS_char.value
ESS_dischar_val = ESS_dischar.value
ESS_val = ESS.value
P_Curt_val = P_Curt.value
P_PV_Output_val = P_PV_Output.value
# print(p.value)
time = list(range(1, 97))

#%% Plot
Fig = plt.figure(figsize=[15,10])
ax1 = plt.subplot(2, 1, 2)
plt.xlabel('Time (h)')
plt.ylabel('OutPut (kW)')
plt.suptitle("PV Curtailment Scheduling", fontsize=25)
plt.plot(time, P_L, 'g-', label = 'Load')
plt.plot(time, P_DG_val, 'r-', label = 'Diesel')
plt.plot(time, ESS_val, 'b-', label = 'ESS')
plt.plot(time, P_PV, 'y-', label = 'PV')
plt.plot(time, P_Curt_val, 'm-', label = 'Curtailment')
plt.plot(time, P_PV_Output_val, 'c-', label = 'PV Output')
plt.legend(loc=2)

plt.subplot(2, 1, 1, sharex=ax1)
plt.ylabel('SOC(%)')
plt.plot(SOC_val, 'b', label = 'SOC')
plt.xticks(visible=False)
plt.legend(loc=2)
plt.show

Fig = plt.figure(figsize=[15,10])
plt.xlabel('time (h)')
plt.ylabel('Output (MW)')
plt.suptitle("Curtailment ", fontsize=25)
plt.plot(P_PV, color='green', label='PV')
plt.plot(P_Curt.value, color='brown', label='Curtail')
plt.plot(P_PV_Output.value, color='blue', label='PV Output')
plt.legend(loc=2)
plt.show
