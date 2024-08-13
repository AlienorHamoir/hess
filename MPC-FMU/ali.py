# Define global simulation parameters
T = 3*60*60  # Time horizon length: horizon T = 3h to start
t_s = 60*60  # Time step size in seconds: ts = 1h -> for not losing too much on transitions for HESS, ideally ts should be smaller or equal to tau_CLD


# Define global parameters
R = 8.314
F = 1.602176634 *10**(-19)
M = 5000
m = 1
V_tn = 1.482
T_amb = 273.15 + 23.5
K=1.4



# Medium parameters


# Electrolyzer parameters
C_rep_EL = 100  # Example replacement cost
N_H_EL = 50     # Example number of hours
C_OM_EL = 10    # Example O&M cost
alpha_EL = C_rep_EL/N_H_EL + C_OM_EL
C_STB_OFF_EL = 5  # Example standby-off cost
C_OFF_STB_EL = 3  # Example off-standby cost
C_ON_OFF_EL = 4   # Example off-on cost
eta_EL = 0.6
eta_F_EL = 0.95 # Faradaic efficiency
P_EL_STB = 285.7
P_EL_max = 5000
P_EL_min = 0.1 * P_EL_max
r_EL = 20 * t_s # 20 W/s 
T_CLD_EL = 1 # minimum interval duration of OFF-STB transition
T_H2,nom = 273.15 + 40
T_EL,nom = 273.15 + 50
N_EL = 20
R_th = 1/22

# Electrolyzer auxiliaries parameters
e_dry,spec = 1400 * 3600 * t_s? # Ws/kg
COP_cw = 4
h_H2 = cp_H2 * (T_H2_nom - T_amb)
h_O2 = cp_O2 * (T_H2_nom - T_amb)
P_pump,disch = 100
P_pump,circ = 285.7
Q_loss_EL = (T_EL,nom - T_amb)/R_th


# Fuel cell parameters
C_rep_FC = 100  # Example replacement cost for the fuel cell
N_H_FC = 50     # Example number of hours for the fuel cell
C_OM_FC = 10    # Example O&M cost for the fuel cell
alpha_FC = C_rep_FC/N_H_FC + C_OM_FC
C_STB_OFF_FC = 5  # Example standby-off cost for the fuel cell
C_OFF_STB_FC = 3  # Example off-standby cost for the fuel cell
C_ON_OFF_FC = 4   # Example on-off cost for the fuel cell
eta_FC = 0.4
eta_F_FC = 0.95 # Faradaic efficiency
P_FC_STB = 300
P_FC_max = 4300
P_FC_min = 300
r_FC = 5 * t_s # 5 W/s 

T_CLD_FC = 1 # minimum interval duration of OFF-STB transition
N_FC = 36
T_Air_in = 273.15 + 23.5
T_H2_in = 273.15 + 40
T_FC,nom = 273.15 + 72
lambda_H2 = 1.5
lambda_O2 = 2.05

# Fuel cell auxiliaries parameters
COP_cw = 4
h_H2_in = cp_H2 * (T_H2_nom - T_amb) # Delta T or just T_in?
h_H2_out = cp_H2 * (T_FC_nom - T_amb)
h_air_in = cp_air * (T_amb - T_amb) 
h_air_out = cp_air * (T_FC_nom - T_amb)
Q_loss_FC = (T_FC_nom - T_amb)/R_th
Q_loss_FC = (T_FC,nom - T_amb)/R_th
eta_comp_air = 0.72
p_comp_in = 1.1325 *10**5 # 1 atm = 1.1325 bar = 1 * 10^5 Pa
p_comp_out = 3 *10**5 # pressure at cathode = 3 bar


# Hydrogen storage parameters
P_hess_nom = 5000 # W
E_hess_nom = 402 # Ws or J
n_H2,max = 21510 #mol/s
SOE_hess_min = 0.2 * E_hess_nom # Ws or J
SOE_hess_max = 0.9 * E_hess_nom # Ws or J
eta_comp = 0.72
p_max = 350*10**5 # 350 bar = 350 * 10^5 Pa
p_comp_in = 1.1325 *10**5 # 1 atm = 1.1325 bar = 1 * 10^5 Pa

LOH_min = 0.1
LOH_max = 0.99

# Battery parameters
E_battery_nom = 5000
V_nom = 230
C_nom = 5000
SOC_battery_min = 0.1
SOC_battery_max = 0.9
SOE_battery_min = 0.2 * E_battery_nom
SOE_battery_max = 0.9 * E_battery_nom
C_battery_inv = 1000  # Example investment cost for the battery
N_battery_cycles = 5000  # Example number of battery cycles

P_battery_nom = 5000
c_battery = C_battery_inv / (N_battery_cycles * P_battery_nom)
eta_batt_ch = 0.98
eta_batt_disch = 0.98
SOE_min = SOE_battery_min
SOE_max = SOE_battery_max

# Grid arameters
P_grid_max = 5000



def optimizer(T, t_s, c_el, pv, load, SOE_initial, LOH_initial): 
    # adapt inputs: T (horizon of simulation), ts (sampling time), SOE_battery, LOH_hess, pv, load, c_el (electricity prices)
    # states: SOE_battery, LOH_hess
    # inputs to FMU: P_battery, P_FC_sys (start with P_FC), P_El_sys (start with P_El), state_FC, state_EL

    # ITERATION 1: no auxiliaries, no thermal model, no transitions
    # inputs variables in MPC: P_battery_disch, P_battery_ch, P_grid, P_FC, P_FC_sys, P_FC_aux, delta_ON_FC, delta_OFF_FC, delta_STB_FC, P_EL, P_EL_sys, P_EL_aux, delta_ON_EL, delta_OFF_EL, delta_STB_EL
    
    
    # Access global parameters
    global P_grid_max, P_battery_nom, eta_batt_ch, eta_batt_disch, SOE_min, SOE_max
    global LOH_min, LOH_max, eta_El, eta_FC, P_EL_min, P_EL_max, P_EL_STB, P_FC_min, P_FC_max, P_FC_STB, r_EL, r_FC
    global alpha_EL, alpha_FC, c_battery

    N = T / t_s
    
    # Initialize
    import casadi
    opti = casadi.Opti()

    # -----------------------------
    # Variables and solver
    # -----------------------------

    SOE = opti.variable(1,N+1)  # state
    LOH = opti.variable(1,N+1)  # state
    P_batt_disch = opti.variable(1,N)       # input
    P_batt_ch = opti.variable(1,N)   # input
    P_grid = opti.variable(1,N) # input 
    P_FC_sys = opti.variable(1,N) # input
    P_FC_aux = opti.variable(1,N) # input
    P_FC = opti.variable(1,N) # input
    delta_ON_FC = opti.variable(1,N) # input
    delta_STB_FC = opti.variable(1,N) # input
    delta_OFF_FC = opti.variable(1,N) # input
    P_EL_sys = opti.variable(1,N) # input
    P_EL_aux = opti.variable(1,N) # input
    P_EL = opti.variable(1,N) # input
    delta_ON_EL = opti.variable(1,N) # input
    delta_STB_EL = opti.variable(1,N) # input
    delta_OFF_EL = opti.variable(1,N) # input

    # delta_HP is a discrete variable (binary)
    discrete_var = [0]*(N+1) + [0]*(N+1) + [0]*N + [0]*N + [0]*N + [0]*N + [0]*N + [0]*N + [1]*N + [1]*N + [1]*N + [0]*N + [0]*N + [0]*N + [1]*N + [1]*N +  [1]*N 

    # Solver
    opti.solver('bonmin', {'discrete': discrete_var, 'bonmin.tol': 1e-4, 'bonmin.print_level': 0, 'print_time': 0})

    # -----------------------------
    # Constraints
    # -----------------------------

    # Initial storage level
    opti.subject_to(SOE[0] == SOE_initial)
    opti.subject_to(LOH[0] == LOH_initial)

    # Constraints at every time step
    for t in range(N+1):

        if t < N:
            # Bounds on grid
            opti.subject_to(pv[t] + P_grid[t] + P_batt_disch[t] + P_FC_sys[t] == load[t] + P_batt_ch[t] + P_EL_sys[t])
            opti.subject_to(P_grid[t]>=0)
            opti.subject_to(P_grid[t]<=P_grid_max)

            # Bounds on battery and SOE
            opti.subject_to(P_batt_disch[t]>=0)
            opti.subject_to(P_batt_disch[t]<=0.9 * P_battery_nom)
            opti.subject_to(P_batt_ch[t]>=0)
            opti.subject_to(P_batt_ch[t]<=0.9 * P_battery_nom)
            
            opti.subject_to(SOE[t+1] == SOE[t] + (P_batt_ch[t] * eta_batt_ch - P_batt_disch[t]/eta_batt_disch) * t_s
            opti.subject_to(SOE[t] >= SOE_min)
            opti.subject_to(SOE[t] <= SOE_max)
            
            # Bounds on hydrogen storage 
            opti.subject_to(LOH[t+1] == LOH[t] + (P_EL[t] * eta_EL - P_FC[t]/eta_FC) * t_s /C_max
            opti.subject_to(SOE[t] >= LOH_min)
            opti.subject_to(SOE[t] <= LOH_max)                

            # Constraints on electrolyzer
            # Bounds when ON                
            opti.subject_to(P_EL[t] >= delta_ON_EL[t] * P_EL_min)
            opti.subject_to(P_EL[t] <= delta_ON_EL[t] * P_EL_max)
            # Bounds when STB
            opti.subject_to(P_EL_aux[t] == delta_STB_EL[t] * P_EL_STB)
            # Ramping constraints                
            opti.subject_to(P_EL[t+1] - P_EL[t] <= r_EL)
            opti.subject_to(P_EL[t] - P_EL[t+1] <= r_EL)
                            
            # Constraints on fuel cell
            # Bounds when ON
            opti.subject_to(P_FC[t] >= delta_ON_FC[t] * P_FC_min)
            opti.subject_to(P_FC[t] <= delta_ON_FC[t] * P_FC_max)
            # Bounds when STB
            opti.subject_to(P_FC_aux[t] == delta_STB_FC[t] * P_FC_STB)
            # Ramping constraints                
            opti.subject_to(P_FC[t+1] - P_FC[t] <= r_FC)
            opti.subject_to(P_FC[t] - P_FC[t+1] <= r_FC)
                            
            # States                
            opti.subject_to(delta_ON_EL[t] + delta_STB_EL[t] + delta_OFF_EL[t] == 1)
            opti.subject_to(delta_ON_FC[t] + delta_STB_FC[t] + delta_OFF_FC[t] == 1)

            # System                
            opti.subject_to(P_EL_sys[t] ==  P_EL[t] + P_EL_aux[t] )
            opti.subject_to(P_FC_sys[t] ==  P_FC[t] - P_FC_aux[t] )
                                          
                            
            # Bilinear to linear


    # -----------------------------
    # Objective
    # -----------------------------

    obj = sum(alpha_EL * delta_ON_EL[t] + alpha_FC * delta_ON_FC[t] + c_battery * abs(P_batt_ch[t] - P_batt_disch[t]) + c_el[t]*P_grid[t]*t_s for t in range(N))
    opti.minimize(obj)

    # -----------------------------
    # Solve and get optimal values
    # -----------------------------

    sol = opti.solve()
    
    SOE_opt = sol.value(SOE)
    LOH_opt = sol.value(LOH)
    P_batt_disch_opt = sol.value(P_batt_disch)
    P_batt_ch_opt = sol.value(P_batt_ch)
    P_grid_opt = sol.value(P_grid)
    P_FC_sys_opt = sol.value(P_FC_sys)
    P_FC_aux_opt = sol.value(P_FC_aux)
    P_FC_opt = sol.value(P_FC)  
    delta_ON_FC_opt = sol.value(delta_ON_FC)
    delta_STB_FC_opt = sol.value(delta_STB_FC)
    delta_OFF_FC_opt = sol.value(delta_OFF_FC)
    P_EL_sys_opt = sol.value(P_EL_sys)
    P_EL_aux_opt = sol.value(P_EL_aux)
    P_EL_opt = sol.value(P_EL)
    delta_ON_EL_opt = sol.value(delta_ON_EL)
    delta_STB_EL_opt = sol.value(delta_STB_EL)
    delta_OFF_EL_opt = sol.value(delta_OFF_EL)

            
    obj_opt = round(sol.value(obj)/100,2)

    

    return SOE_opt, LOH_opt, P_batt_disch_opt, P_batt_ch_opt, P_grid_opt, delta_ON_EL_opt, delta_STB_EL_opt, delta_OFF_EL_opt, P_EL_sys_opt, P_EL_aux_opt, P_EL_opt, delta_ON_FC_opt, delta_STB_FC_opt, delta_OFF_FC_opt, P_FC_sys_opt, P_FC_aux_opt, P_FC_opt, obj_opt




import numpy as np
import matplotlib.pyplot as plt



c_el = np.array([100*10**(-6)/3600, 120*10**(-6)/3600, 110*10**(-6)/3600])  # Electricity prices
pv = np.array([3000, 3500, 4000])  # PV generation
load = np.array([4230, 4500, 4230])  # Load
SOE_initial = 2500  # Initial state of energy of the battery
LOH_initial = 500  # Initial level of hydrogen

# Call the optimizer function
results = optimizer(T, t_s, c_el, pv, load, SOE_initial, LOH_initial)

# Unpack the results
(SOE_opt, LOH_opt, P_batt_disch_opt, P_batt_ch_opt, P_grid_opt, delta_ON_EL_opt, delta_STB_EL_opt, delta_OFF_EL_opt, 
 P_EL_sys_opt, P_EL_aux_opt, P_EL_opt, delta_ON_FC_opt, delta_STB_FC_opt, delta_OFF_FC_opt, 
 P_FC_sys_opt, P_FC_aux_opt, P_FC_opt, obj_opt) = results

# Plot the results
time = np.arange(N+1)  # Time vector for plotting

plt.figure(figsize=(12, 8))

# Plot SOE
plt.subplot(2, 2, 1)
plt.plot(time, SOE_opt, marker='o')
plt.title('State of Energy (SOE) of the Battery')
plt.xlabel('Time [hours]')
plt.ylabel('SOE [Wh]')
plt.grid(True)

# Plot LOH
plt.subplot(2, 2, 2)
plt.plot(time, LOH_opt, marker='o')
plt.title('Level of Hydrogen (LOH)')
plt.xlabel('Time [hours]')
plt.ylabel('LOH [Wh]')
plt.grid(True)

# Plot P_grid
plt.subplot(2, 2, 3)
plt.plot(time[:-1], P_grid_opt, marker='o')
plt.title('Power Grid Usage')
plt.xlabel('Time [hours]')
plt.ylabel('Power [W]')
plt.grid(True)

# Plot P_EL and P_FC
plt.subplot(2, 2, 4)
plt.plot(time[:-1], P_EL_opt, marker='o', label='P_EL')
plt.plot(time[:-1], P_FC_opt, marker='o', label='P_FC')
plt.title('Electrolyzer and Fuel Cell Power')
plt.xlabel('Time [hours]')
plt.ylabel('Power [W]')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
