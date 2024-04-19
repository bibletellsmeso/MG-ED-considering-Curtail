# library import
import cvxpy as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy

#%%--------------------------Hawaai Problem with exact parameters-------------#
# Generate data
P_PV = np.loadtxt("PV_worst_case.txt")                                     # 6/7 Solar data Input
P_L = np.loadtxt("Load_worst_case.txt")                                    # 6/7 Load data Input

T = len(P_L)                                                                   # Total time slots as 1 slot is 15 min so in total it will be 24*4
t = 15/60                                                                      # Delta t is 1/4 hr

# DG data
DG_fuel = 1191.84
DG_a = 31.324
DG_b = 8.786
P_rate = 80
DG_OM = 1.91
ES_OM = 0.955
PV_OM = 1.91
ES_loss = 46.585
PV_penalty = 0.0
DG_R_up = 35
DG_R_down = 35

P_DGmin = 200                                                              # Diesel is always on, minimal power limit 30% of maximum
P_DGmax = 750
i = 0
piece = 10
#%%--------------------------Battery parameters-------------------------------#
ESS_max = 250
ESS_CAP = 750
C_eff = 93/100
D_eff = 100/93
#%%------------------------Variables intialization----------------------------#
# Unit
P_DG = cp.Variable(T)                                                          # DG variable

P_Curt = cp.Variable(T)                                                        # PV Curtailment
P_PV_Output = cp.Variable(T)

curt = cp.Variable(T)                                                        # 1 Piece

# ESS 
ESS_char = cp.Variable(T)
ESS_dischar = cp.Variable(T)
ESS = cp.Variable(T)

SOC = cp.Variable(T)
#%%----------------------------constraints------------------------------------#
constraints = []
for i in range(T):
    temp = []
    # Load balance
    temp.append(P_DG[i] - ESS_char[i] + ESS_dischar[i] + P_PV[i] - P_Curt[i] == P_L[i])

    # ESS
    temp.append(ESS_char[i] <= ESS_max)
    temp.append(ESS_char[i] >= 0)
    temp.append(ESS_dischar[i] <= ESS_max)
    temp.append(ESS_dischar[i] >= 0)
    temp.append(ESS[i] == ESS_char[i] - ESS_dischar[i])

    # SOC
    temp.append(SOC[i] >= 20)
    temp.append(SOC[i] <= 80)
    temp.append(SOC[i] == SOC[i-1] + ((ESS_char[i] * 100 * C_eff * t) / (ESS_CAP)) - ((ESS_dischar[i] * 100 * D_eff * t) / (ESS_CAP))
                if i >= 1 else SOC[i] == 50)
    temp.append(SOC[95] == 50)

    # Generation capacity
    temp.append(P_DG[i] >= P_DGmin)
    temp.append(P_DG[i] <= P_DGmax)

for i in range(T-1):
    temp.append(P_DG[i+1] - P_DG[i] <= DG_R_up)
    temp.append(P_DG[i] - P_DG[i+1] <= DG_R_down)

for i in range(T):
    # PV Curtailment
    temp.append(P_Curt[i] >= 0)
    temp.append(P_Curt[i] <= P_PV[i])

    # sum every time
    constraints += temp
#%%---------------------------------objective function------------------------#
obj = 0
for i in range(T):
    obj += t * (PV_OM * P_PV[i] + PV_penalty * P_Curt[i] + DG_fuel * (DG_a * P_DG[i] + DG_b * P_rate) + DG_OM * P_DG[i] + (ES_OM + ES_loss) * (ESS_char[i] + ESS_dischar[i]))
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
print(p.value)

#%% Plot
Fig = plt.figure(figsize=[15,10])
ax1 = plt.subplot(2, 1, 2)
plt.xlabel('Time (h)')
plt.ylabel('OutPut (kW)')
plt.suptitle("PV Curtailment Scheduling", fontsize=25)
plt.plot(P_L, 'g-', label = 'Load')
plt.plot(P_DG, 'r-', label = 'Diesel')
plt.plot(ESS, 'b-', label = 'ESS')
plt.plot(P_PV, 'y-', label = 'PV')
plt.plot(P_Curt, 'm-', label = 'Curtailment')
plt.plot(P_PV_Output, 'c-', label = 'PV Output')
plt.legend(loc=2)
plt.show()

plt.subplot(2, 1, 1, sharex=ax1)
plt.ylabel('SOC(%)')
plt.plot(SOC_val, 'b', label = 'SOC')
plt.xticks(visible=False)
plt.legend(loc=2)
plt.show()

Fig = plt.figure(figsize=[15,10])
plt.xlabel('time (h)')
plt.ylabel('Output (MW)')
plt.suptitle("Curtailment ", fontsize=25)
plt.plot(P_PV, color='green', label='PV')
plt.plot(P_Curt.value, color='brown', label='Curtail')
plt.plot(P_PV_Output.value, color='blue', label='PV Output')
plt.legend(loc=2)
plt.show()
