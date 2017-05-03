import math
#1. gas constants and ion concentrations
R = 8314.0                                                      #gas constant                          
T = 306.15                                                      #temperature
F = 96487.0                                                     #farady constant
RTF_con = R * T / F                                             #define as RTF constant   
FRT_con = F/(R*T)
 
Na_o = outside_sodium_concentration = 140.0                  #outside_sodium_concentration      
Na_i = inside_sodium_concentration  =18.0                      #inside_sodium_concentration

K_o = outside_potassium_concentration = 5.4                  #outside_potassium_concentration
K_i = inside_potassium_concentration = 145.0                 #inside_potassium_concentration

Ca_o = outside_calcium_concentration = 1.8                   #outside_calcium_concentration
Ca_i = inside_calcium_concentration = 0.00012          #inside_calcium_concentration
#########################____________________finish

#2. calculate Equilibrium
def E(outside,inside):
    '''
    This block is to calculate the Ek of ions.
    Only for mono-valent ions.
    It takes in 2 arguments.
        (1).The first  arg is the EXtracellular concentration of the ion.
        (2).The second arg is the INtracellular concentration of the ion
    For example, in order to calculate the equilibrium potential of sodium channel, try:
        E(137,8)
    '''
    import math
    R = 8314.0
    T = 306.15
    F = 96487.0
    rtf_con = R*T/F
    rate_con = 1.0 * outside/inside
    Ek = rtf_con*math.log(rate_con)
    return Ek
#########################____________________finish

#3.stimulating parameters.
'''
this block try to build the Iapp
I_duration = deltaT
I_amplitude= I[i]            #lower_i is not current, but counter
'''
deltaT=0.01                  #deltaT=0.01 S
t = []                       #set t
for ZZ in range(1000000):     #t = sep(0,1000,0.01)
    t.append(ZZ/100.0)
I = []
for i in range(1000000):     #set I as t
    I.append(0)
for i in range(500):         #blank the beginning interval
    I[i] = 0
for i in range(50500,50700): #begin to stimulate
    I[i] = 20

for i in range(125500,125700): 
    I[i] = 80
#########################____________________finish

#4. Equilibrium Value and G_ion(conductance)
#4a: E_Na_complex   
gbar_Na = 16
E_Na    = E(Na_o,Na_i)

#4b: E_Ca_complex
gbar_Ca = 0.9
E_si    = 7.7 - 13.0827 * math.log(Ca_i)

#4c: E_K_complex
gbar_K  = 0.282  * math.sqrt(K_o/5.4)
PR_NaK = 0.01833
E_K_numerator   = K_o + PR_NaK * Na_o
E_K_denominator = K_i + PR_NaK * Na_i
E_K = RTF_con * math.log(E_K_numerator/E_K_denominator)

#4d: E_K1_complex
gbar_K1 = 0.75 * math.sqrt(K_o/5.4)



gbar_b  = 0.03921
gbar_Kp = 0.0183 

#4e: E_b_complex
E_K1    = E(K_o,K_i)
E_b     = -59.87
E_Kp    = E(K_o,K_i)

#5. parameters of activations and inactivations
#5a. I_Na
def alpha_m(V):
    import math
    numerator   = 0.32 * (V + 47.13)
    denominator =   1  - math.e**(-0.1 * (V + 47.13) )
    return numerator / denominator
def beta_m(V):
    import math
    return 0.08 * math.e**( -V / 11.0)
def alpha_h(V):
    import math
    if V >= -40:
        return 0
    elif V < -40:
        return 0.135 * math.e**((80 + V) / -6.8)
def beta_h(V):
    import math
    if V >= -40:
        numerator_1   = 1.0 
        denominator_1 = 0.13 * (1 + math.e**((V + 10.66)/-11.1))
        return 1 / denominator_1
    elif V < -40:
        component_1 = 3.56 * math.e**(0.079 * V) 
        component_2 = 3.1 * 100000 * math.e**(0.35 * V)
        return component_1 + component_2
def alpha_j(V):
    import math
    if V >= -40:
        return 0
    elif V < -40:
        component_1 = -1.2714 * 100000  * math.e**(0.2444 * V)
        component_2 = 3.474   * 0.00001 * math.e**( - 0.04391 * V)
        denominator =     1   + math.e **(0.311 * (V + 79.23))
        result = (component_1 - component_2) * (V + 37.78)/ denominator
        return result
def beta_j(V):
    import math
    if V >= -40:
        numerator   = 0.3 * math.e**(-2.535 * 0.0000001 * V)
        denominator =   1 + math.e**(-0.1 * (V + 32.0))
        return numerator/denominator
    elif V < -40:
        numerator   = 0.1212 * math.e**(-0.01052 *  V)
        denominator =  1     + math.e**(-0.1378  * (V + 40.14))
        return numerator / denominator
#########################____________________finish

#5b. I_Ca
K_m_Ca_1 = 0.0006    #This should be carefully examined
    
P_Ca_Ca = 0.00054
Gama_Cai = 1.0
Gama_Cao = 0.341

P_Ca_Na = 0.000000675 
Gama_Nai = 0.75
Gama_Nao = 0.75

P_Ca_K  = 0.000000193
Gama_Ki = 0.75
Gama_Ko = 0.75

def inf_d(V):
    import math
    numerator   =    1
    denominator =    1  + math.e**( -(V +10 )/6.24)
    return numerator / denominator

def tau_d(V):
    import math
    numerator   =    1  - math.e**( -(V +10 )/6.24)
    component1  =    1  + math.e**( -(V +10 )/6.24)
    component2  = 0.035 * (V + 10)
    return numerator / ( component1 * component2)

def inf_f(V):
    import math
    component1  =    1  + math.e**( (V  + 35.06 )/8.6 )
    component2  =    1  + math.e**( (50 -   V   )/20.0)
    return (1.0 /  component1) + (0.6/ component2)

def tau_f(V):
    import math
    denominator = 0.0197 * math.e**(-(0.0337 * (V + 10))**2) + 0.02
    return 1 / denominator    

def f_Ca(V):
    import math
    denominator = 1 + (Ca_i / K_m_Ca_1)**2 
    return 1.0 / denominator

def I_general(V,P,z,Gama_i,Gama_o,S_i,S_o):
    #the parameter is V,P_I,z, Gama_in, Gama_out
    R = 8314.0                                                      #gas constant                          
    T = 306.15                                                      #temperature
    F = 96487.0                                                     #farady constant
    RTF_con = R * T / F                                             #define as RTF constant   
    FRT_con = F/(R*T)    
    import math
    component1 = FRT_con * V * F
    interVariable = math.e**(z * V * FRT_con)
    numberator  = Gama_i * S_i * interVariable- Gama_o * S_o
    denominator = interVariable - 1
    return P * z**2 * component1 * numberator / denominator
#########################____________________finish

#5c. I_K
#Warning! This is changed compared to the first model.
def alpha_X(V):
    import math
    numerator   = 0.0000719 * (V + 30)
    denominator = 1  - math.e**(-0.148 * (V + 30))
    return numerator / denominator

def beta_X(V): 
    import math
    numerator   = 0.000131  * (V + 30)
    denominator = -1 + math.e**(0.0687 * (V + 30))
    return numerator / denominator

def Xi(V):
    import math
    numerator   = 1.0
    denominator = 1 +  math.e**((V - 56.26)/32.1)
    return numerator / denominator
#########################____________________finish

#5d. I_K1
def alpha_K1(V):
    import math
    denominator =1 + math.e**(0.2385 * (V - E_K1 - 59.215)) 
    return 1.02 / denominator

def beta_K1(V):
    import math
    component = 0.49124 * math.e**(0.08032 * (V-E_K1+ 5.476))
    enumerator = math.e**(0.06175 * (V- E_K1 - 594.31))
    denominator =1+ math.e**(-0.5143 * (V - E_K1 + 4.753))
    return (component + enumerator) / denominator
#########################____________________finish

#5e. I_Kp
def Kp(V):
    import math
    denominator = 1+ math.e**((7.488 - V) / 5.98)
    return 1.0 / denominator
#########################____________________finish

#5f. I_EX_Na_Ca
def I_EX_Na_Ca(V):
    K_Na_Ca =  2000
    K_m_Na  =  87.5
    K_m_Ca_2 = 1.38
    K_sat    = 0.1
    yita     = 0.35
    component1 = (K_m_Na**3 + Na_o**3) * (K_m_Ca_2 + Ca_o) * (1 + K_sat * math.e**((yita-1) * V * FRT_con))
    component2 = math.e**(yita * V * FRT_con) * Na_i**3 * Ca_o
    component3 = math.e**((yita-1) * V * FRT_con)* Na_o**3 * Ca_i
    return K_Na_Ca / component1 * component2 - component3

#5g. I_pump_Na_K
def I_pump_Na_K(V):
    I_Na_K = 1.5
    K_m_Na_i = 10.0
    K_m_K_o = 1.5
    theta = 1.0 * (math.e**( Na_o/ 67.3) - 1)/ 7  
    component1 = 1 + 0.1245 * math.e**(-0.1 * V * FRT_con) + 0.0365 * theta * math.e**( -V * FRT_con)
    component2 = K_o / ((1 + (K_m_Nai/Na_1)**1.5) * (K_o + K_m_K_o))
    return I_Na_K * component2 / component1

#5h. I_ns_Ca
def I_ns_Ca(V):
    pass


def Ib(V):
    return gbar_b * (V - E_b)
'''
testlist = range(-100,60)

import testplot

'''

