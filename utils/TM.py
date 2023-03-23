import math
#import numpy as np
#import CoolProp.CoolProp as CP


#genral
g = 9.81 #[m/s^2]
delta_p_ps = 47000000.0 #[Pa] allowable pressure rise over a single pump stage (Humble p.254)
N_s = 3.0 #[rad/s] stage specific speed (Humble p.254)
psi = 0.55 # pump head coefficient (Humble p.254)
phi = 0.1 # inducer inlet flow coefficient (Humble, p. 254)
L = 0.3 # inducer inlet hub to tip diameter ratio (Humble, p. 254)



#inputs######################################################

#CEA
p_pb = 30000000  #[Pa]
of_ratio_pb= 0.38 # oxidizer fuel ratio
T_i_t2 = 1169 #[K] (nice) turbine inlet temperature
cp_t2= 9211 #[J/kg*K] constant pressure specific heat
gamma_t2= 1.1361 #C

#pressures
p_cc = 22000000 #[Pa]
p_tank_fuel = 600000 #[Pa]
p_tank_lox = 600000 #[Pa]

#losses
delta_p_injector = 0.1*p_cc #[Pa]
delta_p_cool = 3300000  #[Pa]
delta_p_lines = 500000  #[Pa]

#mass flows
m_dot_fuel = 26.315 #[kg/s] massflow fuel tank
m_dot_ox_tank = 73.685 #[kg/s] massflow tank
m_dot_ox_pb = of_ratio_pb*m_dot_fuel #[kg/s]


#Fuel pump (p3)
T_p3 = 100 #[K] tank temperature                                                               
p_i_p3 = p_tank_fuel #[Pa] pump inlet pressure, no losses in pipes
delta_p_p_p3 = p_pb+delta_p_cool+delta_p_lines-p_i_p3 #[Pa] pump pressure rise                   
p_v_p3 = 116000.0 #[Pa]  vapor pressure                                             ?not in CoolProp

# LOx pump (p4)
T_p4= 90 #[K]
p_i_p4 = p_tank_lox #[Pa] pump inlet pressure 
delta_p_p_p4 = p_cc+delta_p_lines+delta_p_injector-p_i_p4 #[Pa] pump pressure rise                  ?
#gamma_p4 = 1.2 #Annahme
p_v_p4 = 146000.0 #[Pa]  vapor pressure                                             ? not in CoolProp

#LOx boostpump (p2)
delta_p_p_p2 = p_pb+delta_p_lines - p_i_p4 + delta_p_p_p4 #[Pa] pump pressure rise                                   ?
T_p2= 200 #[K]                                                                              !                                                                 
m_dot_p2= m_dot_ox_pb #[kg/s] massflow tank
p_i_p2 = p_i_p4+delta_p_p_p4 #[Pa] pump inlet pressure                                         
p_v_p2 = 146000.0 #[Pa]  vapor pressure                                             ?not in CoolProp
                    
#Turbine (t2)
eta_t2 = 0.8 #Humble p.213
p_i_t2 = p_pb #[Pa] inlet pressure
#P_trat_t2 = 1.5 #turbine pressure ratio between 1,5 an 2,5                         




#calculations##########################################################

# Fuel pump (p3)

roh_p3 = 424 #CP.PropsSI(D,T, T_p3, P, p_i_p3, "methane") #[kg/m^3]

Q_p3 = m_dot_fuel/roh_p3 #[m^3/s] volume flow rate
H_p_p3 = delta_p_p_p3/(g*roh_p3) #[m] pump head rise
NPSH_p3 = (p_i_p3-p_v_p3)/(g*roh_p3) #[m] net positive suction head
n_p3 = math.ceil(delta_p_p_p3/delta_p_ps) # number of stages 
N_r_p3 = (N_s*(H_p_p3/n_p3)**0.75 / math.sqrt(Q_p3)) # [rad/s] pump rotational speed
N_p3 =  (30.0*N_r_p3)/math.pi #[RPM] pump rotational speed
u_t_p3 = math.sqrt((g*H_p_p3)/(n_p3*psi)) #[m/s] pump impeller tip speed
D_2t_p3 = (2*u_t_p3)/N_r_p3 #[m] pum impeller exit tip diameter
D_1t_p3 = (((4/math.pi)*Q_p3)/(phi*N_r_p3*(1-L**2)))**(1/3) #[m] pump impeller inlet tip diameter = inducer exit tip diameter
eta_p_p3 = 0.8 #figure
P_req_p3 = (g*m_dot_fuel*H_p_p3)/eta_p_p3 #[W] power requirement to drive pump


#Lox Pump (p4)
                                            
roh_p4 = 1141.0  # CP.PropsSI("D","T", T_p4, "P", p_i_p4, "oxygen") #[kg/m^3]
m_dot_p4= m_dot_p2 #[kg/s] massflow tank 

Q_p4 = m_dot_p4/roh_p4 #[m^3/s] volume flow rate
H_p_p4 = delta_p_p_p4/(g*roh_p4) #[m] pump head rise
NPSH_p4 = (p_i_p4-p_v_p4)/(g*roh_p4) #[m] net positive suction head
n_p4 = math.ceil(delta_p_p_p4/delta_p_ps) # number of stages,
N_r_p4 = N_r_p3 #(N_s*(H_p_p4/n_p4)**0.75 / math.sqrt(Q_p4)) # [rad/s] pump rotational speed, same shaft no gears
N_p4 =  (30.0*N_r_p4)/math.pi #[RPM] pump rotational speed
u_t_p4 = math.sqrt((g*H_p_p4)/(n_p4*psi)) #[m/s] pump impeller tip speed
D_2t_p4 = (2*u_t_p4)/N_r_p4 #[m] pum impeller exit tip diameter
D_1t_p4 = (((4/math.pi)*Q_p4)/(phi*N_r_p4*(1-L**2)))**(1/3) #[m] pump impeller inlet tip diameter = inducer exit tip diameter
eta_p_p4 = 0.8 #figure
P_req_p4 = (g*m_dot_p4*H_p_p4)/eta_p_p4 #[W] power requirement to drive pump


# LOx boostpump (p2)

roh_p2 = 1141.0 # CP.PropsSI("D","T", T_p2, "P", p_i_p2, "oxygen") #[kg/m^3]

Q_p2 = m_dot_p2/roh_p2 #[m^3/s] volume flow rate
H_p_p2 = delta_p_p_p2/(g*roh_p2) #[m] pump head rise
NPSH_p2 = (p_i_p2-p_v_p2)/(g*roh_p2) #[m] net positive suction head
n_p2 = math.ceil(delta_p_p_p2/delta_p_ps) # number of stages 
N_r_p2 = N_r_p3 #(N_s*(H_p_p2/n_p2)**0.75 / math.sqrt(Q_p2)) # [rad/s] pump rotational speed, same shaft no gears
N_p2 =  (30.0*N_r_p2)/math.pi #[RPM] pump rotational speed
u_t_p2 = math.sqrt((g*H_p_p2)/(n_p2*psi)) #[m/s] pump impeller tip speed
D_2t_p2 = (2*u_t_p2)/N_r_p2 #[m] pum impeller exit tip diameter
D_1t_p2 = (((4/math.pi)*Q_p2)/(phi*N_r_p2*(1-L**2)))**(1/3) #[m] pump impeller inlet tip diameter = inducer exit tip diameter
eta_p_p2 = 0.8                                                                                                                    #figure
P_req_p2 = (g*m_dot_p2*H_p_p2)/eta_p_p2 #[W] power requirement to drive pump


#Turbine (t2)

m_dot_pb = m_dot_fuel + m_dot_ox_pb
P_req_turb = 1.2*(P_req_p2+P_req_p3+P_req_p4)/0.98

P_trat_t2=(1-(P_req_turb/(eta_t2*m_dot_pb*cp_t2*T_i_t2)))**(gamma_t2/(gamma_t2-1)) #Engine balance (Humble p.212)
p_out_t2 = P_trat_t2*p_i_t2

#C_0_t2 = math.sqrt(2*cp_t2*T_i_t2*(1-(1/P_trat_t2))**((gamma_t2-1)/gamma_t2)) #[m/s] isentropic spouting velocity
#u_m_t2 = 425 #[m/s] from diagram maxmimum turbine pitch velocity                                                        !
#vratio_t2 = u_m_t2/C_0_t2 
#D_m_t2 = (2*u_m_t2)/N_r_p3 #mean turbine pitch diameter, without gears


print("Pressure output turbine","%.2f"%(p_out_t2/100000)) #must be bigger than 270
print("Diameter", "%.2f"%D_1t_p2, "%.2f"%D_1t_p3, "%.2f"%D_1t_p4)
print("RPM:","%.0f"%N_p2)
print("Required Power from Turbine:", "%.2f"%P_req_turb)
print("Pressure ratio tubrine:", "%.2f"%P_trat_t2)