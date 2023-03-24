import numpy as np
import CoolProp.CoolProp as CP
import csv
import scipy
import os.path
inc_len = 0.001 #increment inc_length [m]
wt = 0.002 #wall thickness copper liner [m]
depth = 0.001 #depth of a chooling channel [m]
width = 0.0012 #width of a chooling channel [m]
wibc = 0.0005 #width in between channels [m]


#General stuff
lambda_cu = 342 #thermal conductivity copper at 927°C (~1200 K) [W/(m K)], (https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html)
M_CO = 28.01/1000 #[kg/mol]
M_CO2 = 44.01/1000 #[kg/mol]
M_H2O = 18.016/1000 #[kg/mol]
M_H2 = 2.016/1000 #[kg/mol]
M_OH = 17.008/1000 #[kg/mol]
T_c = 100 #Temperature fuel, assuming almost minimally cold liq methane (boiling point ~111K at 1 atm) [K] (https://www.engineeringtoolbox.com/methane-d_1420.html)
cp_m_c = 52.344 #Specific Molar Heat at 100 K, 300 bar [J/(mol K)] (https://webbook.nist.gov/cgi/fluid.cgi?P=30&TLow=100&THigh=150&TInc=25&Digits=5&ID=C74828&Action=Load&Type=IsoBar&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF)
M_CH4 = 16.043 #Molar Mass [g/mol] (https://wtt-pro.nist.gov/wtt-pro/index.html?cmp=methane#methane/A;0,0,507,417;help,about/aa2;50,50,507,417/aa1;475,181,507,417/)
cp_c = cp_m_c*1000 / M_CH4 # [J/(kg K)]
cp_c_w = 4708 #[J/(kg K)] at 1100K (https://www.engineeringtoolbox.com/methane-d_980.html)
p_c = 300 #initial pressure cooling channels [bar]

#From Oscar
x = [] #inc_length along chamber rotational axis[m]
r_nozzle = [] #local radius of nozzle[m]
T_hg =[] #local hot gas temperature [K]
Ma_hg = [] #local hot gas Mach-number 
v_hg = [] #local hot gas velocity [m/s]
dA = [] #local nozzel area between two x increments [m^2]
l = [] #arcinc_length of nozzle surface [m]
token = 1
base = os.path.dirname(__file__)
with open(os.path.join(base, "Files", "flow_data.csv")) as file:
    csvreader = csv.reader(file)

    for row in csvreader:
        if token == 0:
            x.append(float(row[0]))
            r_nozzle.append(float(row[1]))
            T_hg.append(float(row[2]))
            Ma_hg.append(float(row[3]))
            v_hg.append(float(row[4]))
            dA.append(float(row[5]))
            l.append(float(row[6]))
        token = 0

index_t = x.index(0.0)
index_max = len(x)-1
T_t = T_hg[index_t] #temperature at throat [K]
r_t = r_nozzle[index_t] #radius at throat [m]

A_nozzle = 0
for i in range(0, len(dA)):
    A_nozzle += dA[i] #[m^2]
l_arc_nozzel = l[len(l)-1] #inc_length of the nozzle at its surface[m]


#From CEA
rho_cc = 1.4418 #[kg/m^3]
rho_t = 8.9038 #[kg/m^3]
rho_exit = 4.8584 #[kg/m^3]
gamma_cc = 1.1539
gamma_t = 1.1556
gamma_exit = 1.2196
a_cc = 1335.7 #sonic velocity combustion chamber [m/s]
a_t = 1287.8#sonic velocity throat [m/s]
O_F = 2.8 #oxidizer to fuel ratio
cp_hg_t = 4379.4 #cp hot gas at throat [J/(kg*K)]
cp_hg_cc = 4825.7 #cp hot gas in combustion chamber [J/(kg*K)]
cp_hg_exit = 2272.6 #cp hot gas in combustion chamber [J/(kg*K)]
p_cc = 222.91 #Combustion chamber pressure[bar]
p_t = 127.77 #Throat pressure [bar]
M_hg_t = 0.33645*M_CO + 0.19318*M_CO2 + 0.42834*M_H2O + 0.01607*M_H2 + 0.02103*M_OH #Molar mass of hot gas (all mass fractions >1%)
CR_t = 1 #contraction ratio at nozzle
c_star = 1944.1 #[m/s]
p_exit = 0.25 #Nozzle exit pressure [bar]

m_hg = 100 #mass flux hot gas side [kg/s] -> derived from thrust requirement
T_wg_max = 1200 #[K] -> ~150 K from coppers melting point
Ma_t = 1 #Mach number at throat
v_t = Ma_t * a_t #Velocity in throat


#Convert metric to freedom units
rho_t_murica = rho_t *2.2046/(1000000*0.0610) #kg->lb; 1/m^3->1/cm^3->1/in^3 [lb/in^3] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
v_t_murica = v_t *100*0.3937 #m->cm->in [in/s] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
g_murica = 32.2 #earth's gravitational acceleration  [ft/s^2]
M_hg_t_murica = M_hg_t * 2.2046 #[lb/mol] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
T_t_murica = T_t * 9/5 #[°R] (https://www.calculatorsoup.com/calculators/conversions/kelvin-to-rankine.php)
cp_hg_t_murica = cp_hg_t*0.000947817 * 1/(9/5) * 1/2.2046 #J->Btu; °K->°R; kg->lb [Btu/(lb °R)] (https://www.unitconverters.net/energy/joule-to-btu-it.htm)
p_t_murica = 14.503773773 *p_t #[psi] (https://www.unitconverters.net/pressure/bar-to-psi.htm)
r_murica = r_t*100*0.3937 #[in] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
c_star_murica = c_star*3.2808 #[ft/s] (https://www.metric-conversions.org/de/lange/meter-in-fusse.htm)
cp_c_murica = cp_c * 0.000947817 * 1/(9/5) * 1/2.2046 #J->Btu; °K->°R; kg->lb [Btu/(lb °R)] (https://www.unitconverters.net/energy/joule-to-btu-it.htm)




###     START       ####

## GET COOLANT MASS FLUX ##
#Calculations at throat as point of max specific heat flux to get coolant mass flux using boundary condition: Melting point of copper liner T_wg_max

#Heat transfer coefficient hot gas side
Pr_t_apprx = 4*gamma_t / (9*gamma_t - 5) #Prandtl number approximation, Huzel p.86 (eq. 4-15)
sigma = 1/((0.5*(T_wg_max/T_t)*(1+(gamma_t-1)/2 *Ma_t**2)+0.5)**0.68 * (1+(gamma_t-1)/2 * Ma_t**2)**0.12)
D_t = 2*r_t #[m]
D_t_murica = D_t*100*0.3937 #[in]
mu_t_murica = (46.6/(10**10))*M_hg_t_murica**0.5 * T_t_murica**0.6
h_g_t_murica = (0.026/(D_t_murica**0.2)*(mu_t_murica**0.2 * cp_hg_t_murica/(Pr_t_apprx**0.6)) * (p_t_murica*g_murica/c_star_murica)**0.8 * (D_t_murica/r_murica)**0.1)*(CR_t)**0.9 * sigma #Heat transfer coeffiecient hot gas side (Bartz Equation) [Btu/(in^2 s °R)] (Rankine and Fahrenheit have same step size), Huzel Modern Engineering, eq. 4-11, p.85
h_g_t = h_g_t_murica * 20441.748   #Heat transfer coeffiecient hot gas side [W/(m^2 K)] (https://www.translatorscafe.com/unit-converter/de-DE/heat-transfer-coefficient/1-9/watt/meter%C2%B2/K-Btu%20(th)/hour/foot%C2%B2/%C2%B0F/)

#Max specific heat flux:
r = Pr_t_apprx**0.33 #Local recovery factor, Huzel p.85
T_aw_t = T_t* ((1+r*(gamma_t-1)/2 * Ma_t**2)/(1+(gamma_t-1)/2 * Ma_t**2))#adiabatic wall temperature at throat, Huzel p.85 (eq. 4-10-a)
q_max = h_g_t*(T_aw_t - T_wg_max) #[W/m^2]
print(f"{q_max: .2f} W/m^2")


#Coolant side wall tamperature
r2_t = r_t + wt #radius outer side of copper liner at throat [m]
R_th_t_cond = np.log(r2_t/r_t)/(lambda_cu*2*np.pi*inc_len) #Thermal resistance of wall at throat[K/W]; (Péclet equation from WSÜ Arbeitsunterlagen p.17),approximation: cylindrical shell with throat radius
T_wc_t = T_wg_max - q_max*2*np.pi* r_t * inc_len*R_th_t_cond


#Coolant mass flux
mu_c = 174.8/1000000 #[Pa s]
k = 209.1/(10**3) #[W/m K] at 100K and 100 bar closest I could find (https://www.engineeringtoolbox.com/methane-thermal-conductivity-temperature-pressure-d_2021.html)
Pr_c = cp_c*mu_c / k
D_c_single = 2*width*depth / (width+depth) #Hydraulic diameter single cooling channel
n_t = int(2*np.pi*r_t / (width+wibc)) #number of cooling channels at throat
print(f"Number of cooling channels: {n_t}")
D_c = D_c_single *n_t #Hydraulic diameter of all cooling channels
D_c_murica = D_c*100*0.3937#[in]
A_c = width*depth*n_t #crosssection area of all cooling channels [m^2]
A_c_murica = A_c* 10000*0.3937 #[in^2]
A_c_conv = (width*inc_len + inc_len*depth*2)*n_t #surface of convectional heat transfer in cooling channels [m^2]
A_c_conv_murica = A_c_conv * 10000*0.3937 #[in^2]
T_wc_t_murica = T_wc_t * 9/5 #[°R]
T_c_murica = T_c * 9/5 #[°R]
mu_c_murica = 0.4229 * 1/12 * 3600 #dynamic viscosity [lb/(in*s)] (at 100bar, 100 K (closest I could find) https://www.engineeringtoolbox.com/methane-dynamic-kinematic-viscosity-temperature-pressure-d_2068.html)
K = (0.029*cp_c_murica*mu_c_murica**0.2)/(Pr_c**(2/3) * D_c_murica**0.2) * (T_c_murica/T_wc_t_murica)**0.55
Q_max = q_max*r_t *2*np.pi*inc_len #W
Q_max_murica = Q_max * 0.7375621493 *12 #[lb in/s] (https://www.unitconverters.net/power/watt-to-pound-foot-second.htm)
m_c_murica = ((Q_max_murica*A_c_murica**0.8) / (K*A_c_conv_murica*(T_wc_t_murica - T_c_murica)))**(1.25) #[lb/s]
m_c = m_c_murica/2.2046 #[kg/s]
print(f"Purely regenrative cooling mass flow: {m_c: .2f} kg/s")

#Assumption: All methane through cooling channels
m_f = m_hg/(O_F+1) #methane total mass flow[kg/s]
m_f_murica = m_f*2.2046 #[lb/s]
Q_Max_all_mf_murica = (m_f_murica*(K*A_c_conv_murica*(T_wc_t_murica - T_c_murica))**(1.25) / (A_c_murica**0.8))**(1/1.25)
Q_Max_all_mf = Q_Max_all_mf_murica / (0.7375621493 *12) #[W]
q_max_all_mf = Q_Max_all_mf / (r_t *2*np.pi*inc_len) #[W/m^2] Maximum specific heat pick up by regenerative system

#Film cooling calculations
SF = 3 #safety factro for very crude calculation
delta_q = q_max - q_max_all_mf
T_fcool_t = 1000 #[K]  assumption: T_coolant outlet = 550K -> at throat 1000K?
cp_ch4_cc_t = CP.PropsSI('C','P', 128e5,'T',900,'methane') # [J/kgK] 128 bar = pressure at throat (CEA)
m_fcool_specific = delta_q / (cp_ch4_cc_t*(T_wg_max-T_fcool_t)) #[kg/s*m^2]
m_fcool = A_nozzle * m_fcool_specific * SF
print(f"Film cooling mass flow: {m_fcool: .2f} kg/s")



##GETTING COOLANT OUTLET TEMPERATURE##
T_hg_initial = T_hg[index_max]
#Getting the pressure distribution inside of the nozzle/chamber
#Assumption: linear p decrease
p_hg = []
for i in range(0, index_max+1):
    p_hg.append(0)
p_hg[0] = p_cc
p_hg[index_max] = p_exit
p_hg[index_t] = p_t
#for injector until throat:
for i in range(1, index_t):
    lenInjTh = abs(x[0]) #axial length in between throat and injecetor [m], (throat at x =0)
    p_hg[i] = round(p_cc - (p_cc-p_t)/lenInjTh * (x[i]-x[0]), 4) #local hot gas pressure [bar]
#for throat until nozzle exit:
for i in range(index_t+1, index_max):
    lenThExit = x[index_max]  #axial length in between throat and nozzle exit [m], (throat at x =0)
    p_hg[i] = round(p_t - (p_t-p_exit)/lenThExit * (x[i]-x[index_t]), 4) #local hot gas pressure [bar]

#Getting the cp distribution of the nozzle/chamber
#Assumption: linear cp decrease
cp_hg = []
for i in range(0, index_max+1):
    cp_hg.append(0)
cp_hg[0] = cp_hg_cc
cp_hg[index_max] = cp_hg_exit
cp_hg[index_t] = cp_hg_t
#for injector until throat:
for i in range(1, index_t):
    lenInjTh = abs(x[0]) #axial length in between throat and injecetor [m], (throat at x =0)
    cp_hg[i] = round(cp_hg_cc - (cp_hg_cc-cp_hg_t)/lenInjTh * (x[i]-x[0]), 4) #local hot gas pressure [bar]
#for throat until nozzle exit:
for i in range(index_t+1, index_max):
    lenThExit = x[index_max]  #axial length in between throat and nozzle exit [m], (throat at x =0)
    cp_hg[i] = round(cp_hg_t - (cp_hg_t-cp_hg_exit)/lenThExit * (x[i]-x[index_t]), 4) #local hot gas pressure [bar]

#Getting the gamma distribution of the nozzle/chamber
#Assumption: linear gamma increase
gamma_hg = []

for i in range(0, index_max+1):
    gamma_hg.append(0)

gamma_hg[0] = gamma_cc
gamma_hg[index_max] = gamma_exit
gamma_hg[index_t] = gamma_t
#for injector until throat:
for i in range(1, index_t):
    lenInjTh = abs(x[0]) #axial length in between throat and injecetor [m], (throat at x =0)
    gamma_hg[i] = round(gamma_cc - (gamma_cc-gamma_t)/lenInjTh * (x[i]-x[0]), 4) #local hot gas pressure [bar]
#for throat until nozzle exit:
for i in range(index_t+1, index_max):
    lenThExit = x[index_max]  #axial length in between throat and nozzle exit [m], (throat at x =0)
    gamma_hg[i] = round(gamma_t - (gamma_t-gamma_exit)/lenThExit * (x[i]-x[index_t]), 4) #local hot gas pressure [bar]

#Getting the density distribution of the nozzle/chamber
#Assumption: linear density decrease
rho_hg = []

for i in range(0, index_max+1):
    rho_hg.append(0)

rho_hg[0] = rho_cc
rho_hg[index_max] = rho_exit
rho_hg[index_t] = rho_t
#for injector until throat:
for i in range(1, index_t):
    lenInjTh = abs(x[0]) #axial length in between throat and injecetor [m], (throat at x =0)
    rho_hg[i] = round(rho_cc - (rho_cc-rho_t)/lenInjTh * (x[i]-x[0]), 4) #local hot gas pressure [bar]
#for throat until nozzle exit:
for i in range(index_t+1, index_max):
    lenThExit = x[index_max]  #axial length in between throat and nozzle exit [m], (throat at x =0)
    rho_hg[i] = round(rho_t - (rho_t-rho_exit)/lenThExit * (x[i]-x[index_t]), 4) #local hot gas pressure [bar]


p_c_in = p_c
T_c_in = T_c
delta_p = 0


for i in range(index_max, 0, -5):
    #Convert metric to freedom units
    rho_hg_murica = rho_hg[i] *2.2046/(1000000*0.0610) #kg->lb; 1/m^3->1/cm^3->1/in^3 [lb/in^3] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
    v_hg_murica = v_hg[i] *100*0.3937 #m->cm->in [in/s] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
    g_murica = 32.2 #earth's gravitational acceleration  [ft/s^2]
    M_hg_t_murica = M_hg_t * 2.2046 #[lb/mol] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
    T_hg_murica = T_hg[i] * 9/5 #[°R] (https://www.calculatorsoup.com/calculators/conversions/kelvin-to-rankine.php)
    cp_hg_murica = cp_hg[i]*0.000947817 * 1/(9/5) * 1/2.2046 #J->Btu; °K->°R; kg->lb [Btu/(lb °R)] (https://www.unitconverters.net/energy/joule-to-btu-it.htm)
    p_hg_murica = 14.503773773 *p_hg[i] #[psi] (https://www.unitconverters.net/pressure/bar-to-psi.htm)
    r_murica = r_nozzle[i]*100*0.3937 #[in] (https://www.mathsisfun.com/metric-imperial-conversion-charts.html)
    c_star_murica = c_star*3.2808 #[ft/s] (https://www.metric-conversions.org/de/lange/meter-in-fusse.htm)
    cp_c_murica = CP.PropsSI('C','P', p_c_in * 10**5,'T',T_c_in,'methane') * 0.000947817 * 1/(9/5) * 1/2.2046 #J->Btu; °K->°R; kg->lb [Btu/(lb °R)] (https://www.unitconverters.net/energy/joule-to-btu-it.htm)

    #Heat transfer coefficient hot gas side
    Pr_apprx = 4*gamma_hg[i] / (9*gamma_hg[i] - 5) #Prandtl number approximation, Huzel p.86 (eq. 4-15)
    sigma = 1/((0.5*(T_wg_max/T_hg[i])*(1+(gamma_hg[i]-1)/2 *Ma_hg[i]**2)+0.5)**0.68 * (1+(gamma_hg[i]-1)/2 * Ma_hg[i]**2)**0.12)
    D_local = 2*r_nozzle[i] #[m]
    D_local_murica = D_local*100*0.3937 #[in]
    mu_murica = (46.6/(10**10))*M_hg_t_murica**0.5 * T_hg_murica**0.6
    #h_g_murica = (0.026/(D_local_murica**0.2)*(mu_murica**0.2 * cp_hg_murica/(Pr_apprx**0.6)) * (p_hg_murica*g_murica/c_star_murica)**0.8 * (D_local_murica/r_murica)**0.1)*((r_t**2*np.pi)/(r_nozzle[i]**2*np.pi))**0.9 * sigma #Heat transfer coeffiecient hot gas side (Bartz Equation) [Btu/(in^2 s °R)] (Rankine and Fahrenheit have same step size), Huzel Modern Engineering, eq. 4-11, p.85
    h_g_murica = (rho_hg_murica*v_hg_murica)**0.8
    h_g = h_g_murica * 20441.748   #Heat transfer coeffiecient hot gas side [W/(m^2 K)] (https://www.translatorscafe.com/unit-converter/de-DE/heat-transfer-coefficient/1-9/watt/meter%C2%B2/K-Btu%20(th)/hour/foot%C2%B2/%C2%B0F/)

    #Max specific heat flux:
    r = Pr_apprx**0.33 #Local recovery factor, Huzel p.85
    T_aw = T_hg[i]* ((1+r*(gamma_hg[i]-1)/2 * Ma_hg[i]**2)/(1+(gamma_hg[i]-1)/2 * Ma_hg[i]**2))#adiabatic wall temperature at throat, Huzel p.85 (eq. 4-10-a)

    #Coolant side wall tamperature
    r2_nozzle = r_nozzle[i] + wt #radius outer side of copper liner at throat [m]
    R_th_cond = np.log(r2_nozzle/r_nozzle[i])/(lambda_cu*2*np.pi*inc_len) #Thermal resistance of wall at throat[K/W]; (Péclet equation from WSÜ Arbeitsunterlagen p.17),approximation: cylindrical shell with throat radius
    T_wc = T_wg_max - q_max*2*np.pi* r_t * inc_len*R_th_cond

    #Heat transfer coefficient coolant side
    mu_c = CP.PropsSI('V','P', p_c_in * 10**5,'T',T_c_in,'methane') #Dynamic viscoyity [Pa s]
    k = CP.PropsSI('L','P', p_c_in * 10**5,'T',T_c_in,'methane') #Thermal conductivity [W/m K]
    cp_c = CP.PropsSI('C','P', p_c_in * 10**5,'T',T_c_in,'methane') #[J/kg K]
    Pr_c = cp_c*mu_c / k
    D_c_single = 2*width*depth / (width+depth) #Hydraulic diameter single cooling channel
    n_t = int(2*np.pi*r_nozzle[i] / (width+wibc)) #number of cooling channels at throat
    D_c = D_c_single *n_t #Hydraulic diameter of all cooling channels
    D_c_murica = D_c*100*0.3937#[in]
    A_c = width*depth*n_t #crosssection area of all cooling channels [m^2]
    A_c_murica = A_c* 10000*0.3937 #[in^2]
    A_c_conv = (width*inc_len + inc_len*depth*2)*n_t #surface of convectional heat transfer in cooling channels [m^2]
    A_c_conv_murica = A_c_conv * 10000*0.3937 #[in^2]
    T_wc_murica = T_wc * 9/5 #[°R]
    T_c_murica = T_c_in * 9/5 #[°R]
    mu_c_murica = mu_c * 0.05599741 #dynamic viscosity [lb/(in*s)] (https://www.aqua-calc.com/convert/dynamic-viscosity/pascal-second-to-pound-per-inch-second)
    m_f_murica = m_f*2.2046 #[lb/s]

    h_c_murica = 0.029*cp_c*mu_c**0.2 / (Pr_c**(2/3)) * ((m_f_murica/A_c_murica)**0.8 / (D_c_murica**0.2) * (T_c_in/T_wc)**0.55)
    #h_c_murica = 0.7
    h_c = h_c_murica * 20441.748   #Heat transfer coeffiecient coolant side [W/(m^2 K)] (https://www.translatorscafe.com/unit-converter/de-DE/heat-transfer-coefficient/1-9/watt/meter%C2%B2/K-Btu%20(th)/hour/foot%C2%B2/%C2%B0F/)

    #Total heat transfer coefficient
    U = 1/ (1/h_g + r_nozzle[i]/lambda_cu * np.log(r2_nozzle/r_nozzle[i]) + r_nozzle[i]/(r2_nozzle * h_c))
    A = 3.1415926*r_nozzle[i]*2*(x[i] - x[i-1])#[m^2]
    T_hg_avg = abs(T_hg[i] + T_hg[i-1])*0.5
    T_c_out = T_c_in
    while True:
        T_c_avg = abs(T_c_in+T_c_out)*0.5
        Q =  U * A * (T_hg_avg - T_c_avg)
        delta_T = Q/(cp_c*m_f)
        if T_c_out*1.001 > (T_c_in + delta_T) and T_c_out*0.999 < (T_c_in+delta_T):
            break
        T_c_out = T_c_in + delta_T

    #Pressure drop
    Re_c = 4*m_f/(np.pi *mu_c * D_c)
    delta_p = delta_p + abs(0.092*(Re_c**1.8 * (l[i]-l[i-1]) * mu_c**2)/(CP.PropsSI('D','P', p_c_in * 10**5,'T',T_c_in,'methane') * D_c**3))

    T_c_in = T_c_out

print(f"{T_c_out: .2f} K")
#Re_c = 4*m_f/(np.pi *mu_c * D_c)
#delta_p =  0.092*(Re_c**1.8 * l_arc_nozzel * mu_c**2)/(CP.PropsSI('D','P', p_c_in * 10**5,'T',T_c_in,'methane') * D_c**3)
print(delta_p)

        



