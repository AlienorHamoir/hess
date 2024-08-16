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



# Discrete logical states 2

def enforce_state_constraints(k, delta_ON, delta_OFF, delta_STB, z_min_pos, z_max_neg, z_stb_pos, z_stb_neg, z_off_pos, z_off_neg):
    constraints = []

    # 1. (1 - δ^ON_i(k)) + z^>=P^min_i(k) >= 1
    constraints.append((1 - delta_ON) + z_min_pos >= 1)

    # 2. (1 - δ^ON_i(k)) + z^<=P^max_i(k) >= 1
    constraints.append((1 - delta_ON) + z_max_neg >= 1)

    # 3. (1 - δ^STB_i(k)) + z^>=P^STB_i(k) >= 1
    constraints.append((1 - delta_STB) + z_stb_pos >= 1)  

    # 4. (1 - δ^STB_i(k)) + z^<=P^STB_i(k) >= 1
    constraints.append((1 - delta_STB) + z_stb_neg >= 1)  

    # 5. (1 - δ^OFF_i(k)) + z^>=0_i(k) >= 1
    constraints.append((1 - delta_OFF) + z_off_pos >= 1)  

    # 6. (1 - δ^OFF_i(k)) + z^<=0_i(k) >= 1
    constraints.append((1 - delta_OFF) + z_off_neg >= 1) 

    return constraints


# Check/Write? constraints for EL
states_constraints_EL = enforce_state_constraints(
    delta_ON_EL[k], z_min_pos_EL[k], z_max_neg_EL[k], z_stb_pos_EL[k], z_stb_neg_EL[k], z_pos_EL[k], z_neg_EL[k]
)

# Check/Write? constraints for FC
states_constraints_FC = enforce_state_constraints(
    delta_ON_FC[k], z_min_pos_FC[k], z_max_neg_FC[k], z_stb_pos_FC[k], z_stb_neg_FC[k], z_pos_FC[k], z_neg_FC[k]
)

# Output results
print("Constraints for EL:", states_constraints_EL)
print("Constraints for FC:", states_constraints_FC)


# States transitions 1

# Binary decision variables for EL
delta_ON_EL 
delta_OFF_EL
delta_STB_EL

sigma_ON_OFF_EL
sigma_OFF_ON_EL
sigma_STB_ON_EL 
sigma_ON_STB_EL 
sigma_STB_OFF_EL 
sigma_OFF_STB_EL 


# Binary decision variables for FC
delta_ON_FC 
delta_OFF_FC
delta_STB_FC

sigma_ON_OFF_FC
sigma_OFF_ON_FC
sigma_STB_ON_FC 
sigma_ON_STB_FC 
sigma_STB_OFF_FC 
sigma_OFF_STB_FC



def enforce_transition_constraints(
    delta_ON_prev, delta_ON, delta_OFF_prev, delta_OFF, delta_STB_prev, delta_STB,
    sigma_ON_OFF, sigma_OFF_ON, sigma_STB_ON, sigma_ON_STB, sigma_STB_OFF, sigma_OFF_STB
):
    constraints = []

    # Constraints ON to OFF
    constraints.append(-delta_ON_prev + sigma_ON_OFF <= 0)
    constraints.append(-delta_OFF + sigma_ON_OFF <= 0)
    constraints.append(delta_ON_prev + delta_OFF - sigma_ON_OFF <= 1)

    # Constraints OFF TO ON
    constraints.append(-delta_OFF_prev + sigma_OFF_ON <= 0)
    constraints.append(-delta_ON + sigma_OFF_ON <= 0)
    constraints.append(delta_OFF_prev + delta_ON - sigma_OFF_ON <= 1)

    # Constraints STB TO ON
    constraints.append(-delta_STB_prev + sigma_STB_ON <= 0)
    constraints.append(-delta_ON + sigma_STB_ON <= 0)
    constraints.append(delta_STB_prev + delta_ON - sigma_STB_ON <= 1)
    
    # Constraints ON TO STB
    constraints.append(-delta_ON_prev + sigma_ON_STB <= 0)
    constraints.append(-delta_STB + sigma_ON_STB <= 0)
    constraints.append(delta_ON_prev + delta_STB - sigma_ON_STB <= 1)
    
    # Constraints STB TO OFF
    constraints.append(-delta_STB_prev + sigma_STB_OFF <= 0)
    constraints.append(-delta_OFF + sigma_STB_OFF <= 0)
    constraints.append(delta_STB_prev + delta_OFF - sigma_STB_OFF <= 1)

    # Constraints OFF TO STB
    constraints.append(-delta_OFF_prev + sigma_OFF_STB <= 0)
    constraints.append(-delta_STB + sigma_OFF_STB <= 0)
    constraints.append(delta_OFF_prev + delta_STB - sigma_OFF_STB <= 1)

    return constraints


# Check constraints for EL
constraints_EL = enforce_transition_constraints(
    delta_ON_EL[k-1], delta_ON_EL[k],  delta_OFF_EL[k-1],  delta_OFF_EL[k], delta_STB_EL[k-1], delta_STB_EL[k],
    sigma_ON_OFF_EL[k], sigma_OFF_ON_EL[k] = 0, sigma_STB_ON_EL[k], sigma_ON_STB_EL[k], sigma_STB_OFF_EL[k], sigma_OFF_STB_EL[k]
)

# Check constraints for FC
constraints_FC = enforce_transition_constraints(
    delta_ON_FC[k-1], delta_ON_FC[k],  delta_OFF_FC[k-1],  delta_OFF_FC[k], delta_STB_FC[k-1], delta_STB_FC[k],
    sigma_ON_OFF_FC[k], sigma_OFF_ON_FC[k] = 0, sigma_STB_ON_FC[k], sigma_ON_STB_FC[k], sigma_STB_OFF_FC[k], sigma_OFF_STB_FC[k]
)

sigma_OFF_ON_EL = np.zeros(T)
sigma_OFF_ON_FC = np.zeros(T)

# Output results
print("Constraints for EL:", constraints_EL)
print("Constraints for FC:", constraints_FC)

# States transitions 2
# Initialize parameters
T_CLD_EL = 1
T_CLD_FC = 1


def enforce_startup_constraints(sigma_STB_OFF, T_CLD, k):
    # Initialize a list to store constraint satisfaction statuses
    constraints = []

    # Generate tau_CLD as a vector of time steps: tau_CLD = (k+1, k+2, ..., k+T_CLD)
    tau_CLD = np.arange(k + 1, k + T_CLD)

    # First constraint: σ^STB_OFF_i(k) - σ^STB_OFF_i(k-1) <= σ^STB_OFF_i(τ_CLD)
    for tau in tau_CLD:
       # if tau < len(sigma_STB_OFF):
        constraints.append(sigma_STB_OFF[k] - sigma_STB_OFF[k-1] <= sigma_STB_OFF[tau])
    
    # Second constraint: Sum of σ^STB_OFF_i(k-T_CLD) to σ^STB_OFF_i(k) <= T_CLD
    if k >= T_CLD:
        sigma_sum = sum(sigma_STB_OFF[k - T_CLD:k])  # Sum over the period
        constraints.append(sigma_sum <= T_CLD)
    else:
        sigma_sum = sum(sigma_STB_OFF[0:k])  # Sum over the period
        constraints.append(sigma_sum <= T_CLD)
        
    return constraints



# Electrolyzer start-up time
startup_constraints_EL = enforce_startup_constraints(sigma_STB_OFF_EL, T_CLD_EL, k)


# Fuel cell start-up time
startup_constraints_FC = enforce_startup_constraints(sigma_STB_OFF_FC, T_CLD_FC, k)


# Electrolyzer dynamics

# Initialize parameters
GCV_H2 = 
M_H2 =
M_O2 =
M_H2O =
cp_H2 = 
cp_O2 = 
eta_EL = 0.6
eta_F_EL = 
R = 8.314
F = 
N_EL = 20
V_tn = 1.482
R_th = 
T_H2_nom = 273.15 + 40
T_EL_nom = 273.15 + 50
T_amb = 273.15 + 23.5
e_dry_spec = 1400 * 3600 * t_s? # Ws/kg
COP_cw = 4
K=1.4
eta_comp = 0.72
p_max = 350*10**5 # 350 bar = 350 * 10^5 Pa
p_comp_in = 1.1325 *10**5 # 1 atm = 1.1325 bar = 1 * 10^5 Pa
P_EL_STB = 
h_H2 = cp_H2 * (T_H2_nom - T_amb)
h_O2 = cp_O2 * (T_H2_nom - T_amb)
P_pump_disch =
P_pump_circ = 
Q_loss_EL = (T_EL,nom - T_amb)/R_th



# mass 
dot_m_H2_EL[k] = eta_EL * P_EL[k] / GCV_H2
dot_n_H2_EL[k] = dot_m_H2_EL[k] / M_H2
dot_n_O2_EL[k] = 0.5 * dot_n_H2_EL[k] 

# current
I_EL[k] = 2 * F * dot_m_H2_EL[k] / (M_H2 * N_EL * eta_F_EL)

# cooling & temperature
Q_cooling_EL[k] = P_gen[k] + P_pump,disch - H_flow_prod[k] - Q_loss_EL[k]
P_gen[k] = P_EL[k] - N_EL * V_tn * I_EL[k]
H_flow_prod[k] = dot_n_H2_EL[k] * h_H2 + dot_n_O2_EL[k] * h_O2 

# auxiliaries power
P_EL_aux[k] = P_dry[k] + P_pump,circ + P_cooling_EL[k] + P_comp,H2[k]
P_dry[k] = e_dry,spec * dot_m_H2_EL[k] 
P_cooling_EL[k] = Q_cooling_EL[k]/COP_cw
RECAST: P_comp,H2[k] = dot_n_H2_EL[k] / eta_comp * (K/K-1) * R * T_H2,nom * ( (p_max * SOC_hess[k]/p_comp_in)**(K-1/K) - 1 )


# electrolyzer system power
P_EL_sys[k] = z_EL_in[k] + z_EL_aux[k] + P_EL_STB * delta_STB_EL[k]


# Fuel cell dynamics

# Initialize parameters
NCV_H2 = 
M_H2 =
M_O2 =
M_H2O =
cp_H2 = 
cp_air =
cp_H2O = 
Delta_H0 =
Delta_HT = cp_H2O * (T_FC_nom - T_amb)
eta_FC = 0.4
eta_F_FC = 
R = 8.314
F = 
N_FC = 36
V_tn = 1.482
R_th_FC = 
T_Air_in = 273.15 + 23.5
T_H2_in = 273.15 + 40
T_FC,nom = 273.15 + 72
T_amb = 273.15 + 23.5
lambda_H2 = 1.5
lambda_O2 = 2.05
COP_cw = 4
K=1.4
eta_comp_air = 0.72
p_comp_in = 1.1325 *10**5 # 1 atm = 1.1325 bar = 1 * 10^5 Pa
p_comp_out = 3 *10**5 # pressure at cathode = 3 bar
P_FC_STB = 
h_H2_in = cp_H2 * (T_H2,nom - T_amb) # Delta T or just T_in?
h_H2_out = cp_H2 * (T_FC,nom - T_amb)
h_air_in = cp_air * (T_amb - T_amb) 
h_air_out = cp_air * (T_FC,nom - T_amb)
Q_loss_FC = (T_FC,nom - T_amb)/R_th



# mass 
dot_m_H2_react[k] = eta_FC * P_FC[k] / NCV_H2
dot_m_H2_in[k] = dot_m_H2_react[k] * lambda_H2
dot_m_H2_out[k] = dot_m_H2_in[k] - dot_m_H2_react[k]
dot_n_H2_react[k] = dot_m_H2_react[k] / M_H2

dot_n_O2_react[k] = 0.5 * dot_n_H2_react[k] 
dot_n_air_react[k] = dot_n_O2_react[k] / 0.21
dot_m_O2_react[k] = dot_n_O2_react[k] * M_O2
dot_m_O2_in[k] = dot_m_O2_react[k] * lambda_O2
dot_m_air_in[k] = dot_m_O2_in[k] / 0.21
dot_m_O2_out[k] = dot_m_O2_in[k] - dot_m_O2_react[k]
dot_m_air_out[k] = dot_m_O2_out[k]/0.21

dot_m_H2O_gen[k] = dot_n_H2_react[k] * M_H2O

# current
I_FC[k] = 2 * F * eta_F_FC * dot_m_H2_react[k] / (M_H2 * N_FC)

# cooling & temperature
Q_cooling_FC[k] = Q_tot[k] +  Q_sens[k] - P_FC[k] - Q_loss_FC
Q_tot[k] = dot_m_H2_react[k] * Delta_H0
Q_sens[k] = dot_m_H2_in[k] * h_H2_in -  dot_m_H2_out[k] * h_H2_out + dot_m_air_in[k] * h_air_in -  dot_m_air_out[k] * h_air_out 
                - dot_m_H2O_gen[k] * Delta_HT/M_H2O

# auxiliaries power
P_FC_aux[k] = P_cooling_FC[k] + P_comp_air[k] 
P_cooling[k] = Q_cooling_FC[k]/COP_cw
P_comp_air[k] = dot_n_air_EL[k] / eta_comp_air * (K/K-1) * R * T_Air_in * ( (p_comp_out/p_comp_in)**(K-1/K) - 1 )


# electrolyzer system power
P_FC_sys[k] = z_FC_out[k] - z_FC_aux[k] - P_FC_STB * delta_STB_FC[k]


# Electrolyzer and fuel cell operating constraints

# Initialize parameters
P_EL_min = 0.1 * P_EL_max
P_EL_max = 5000
P_FC_min = 300
P_FC_max = 4300
r_EL = 20 * t_s # 20 W/s 
r_FC = 5 * t_s # 5 W/s 


# Min and max power boundaries
P_EL[k] >= P_EL_min
P_EL[k] <= P_EL_max

P_FC[k] >= P_FC_min
P_FC[k] <= P_FC_max


# Ramping constraints
#np.abs((P_EL[k+1]-P_EL[k])*delta_ON_EL[k]) <= r_EL
#np.abs((P_FC[k+1]-P_FC[k])*delta_ON_FC[k]) <= r_FC
np.abs(z_EL_in[k+1]-z_EL_in[k]) <= r_EL
np.abs(z_FC_out[k+1]-z_FC_out[k]) <= r_EL


# Logical powers: recast P_EL(k) * delta_ON_EL(k) and P_FC(k) * delta_ON_FC(k) and replace in the code
z_EL_in[k] >= m * delta_ON_EL[k]
z_EL_in[k] <= M * delta_ON_EL[k]
z_EL_in[k] >= P_EL[k] - M * (1 - delta_ON_EL[k])
z_EL_in[k] <= P_EL[k] + m * (1 - delta_ON_EL[k])


z_FC_out[k] >= m * delta_FC_EL[k]
z_FC_out[k] <= M * delta_FC_EL[k]
z_FC_out[k] >= P_FC[k] - M * (1 - delta_ON_FC[k])
z_FC_out[k] <= P_FC[k] + m * (1 - delta_ON_FC[k])

# Logical powers for auxiliaries: recast P_EL,aux(k) * delta_ON_EL(k) and P_FC,aux(k) * delta_ON_FC(k) and replace in the code
z_EL_aux[k] >= m_aux * delta_ON_EL[k]
z_EL_aux[k] <= M_aux * delta_ON_EL[k]
z_EL_aux[k] >= P_EL[k] - M_aux * (1 - delta_ON_EL[k])
z_EL_aux[k] <= P_EL[k] + m_aux * (1 - delta_ON_EL[k])


z_FC_aux[k] >= m_aux * delta_FC_EL[k]
z_FC_aux[k] <= M_aux * delta_FC_EL[k]
z_FC_aux[k] >= P_FC[k] - M_aux * (1 - delta_ON_FC[k])
z_FC_aux[k] <= P_FC[k] + m_aux * (1 - delta_ON_FC[k])
