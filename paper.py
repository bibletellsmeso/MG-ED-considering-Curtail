import os
# GRB: NumPy와 Scipy 패키지 뿐만 아니라 Gurobi 함수와 클래스 불러오기
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

PV = np.loadtxt("PV_worst_case.txt")                                     # 6/7 Solar data Input
L = np.loadtxt("Load_worst_case.txt")                                    # 6/7 Load data Input
ESS = np.loadtxt("ESS_worst_case.txt")

Sim_time = 96
N_PWL = 10
RTE = 0.93
DG_fuel = 1191.84
DG_a = 31.324
DG_b = 8.786
P_rate = 80
DG_OM = 1.91
ES_OM = 0.955
PV_OM = 1.91
ES_loss = 46.585
PV_penalty = 57.3
DG_R_up = 35
DG_R_down = 35

de = gp.Model("de_problem")

var_pv = PV
var_demand = L
var_ess = ESS
var_dg = de.addMVar(shape=Sim_time, lb=0, ub=750, vtype=GRB.CONTINUOUS)

