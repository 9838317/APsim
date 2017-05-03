#this is to change the dir, and confirm that the program could be working under different OS.
import os
import sys
current_path = os.getcwd()                              #get current dir
parent_path = os.path.dirname(current_path)             #get parent  dir
sys.path.append(current_path)                           
sys.path.append(parent_path)                            #to append the previous

#import constant
from z_pa_and_con import *
import numpy as np

def calculate_alpha(inf, tau):
    return inf / tau
def calculate_beta(inf, tau):
    return (1.0-inf) /tau
def probability_inf(a,b):
    return 1.0 * a / (a + b)


import matplotlib.pyplot as plt 

#first cycle:
#CALCULATE alpha,beta, then m,n
class actionpo(object):
    def __init__(self):
        #Create all the list
        #0.V  
        self.list_V       = []
        #1.I_Na
        self.list_alpha_m = []
        self.list_beta_m  = []
        self.list_alpha_h = []
        self.list_beta_h  = []
        self.list_alpha_j = []
        self.list_beta_j  = []
        self.list_m       = []
        self.list_h       = []
        self.list_j       = []
        self.list_I_Na    = []
        #2.I_Ca
        self.list_alpha_d = []
        self.list_beta_d  = []
        self.list_alpha_f = []
        self.list_beta_f  = []
        self.list_tau_d   = []
        self.list_tau_f   = []

        self.list_inf_d   = []
        self.list_inf_f   = []
        
        self.list_f_Ca    = []
        
        self.list_I_Ca_Ca = []
        self.list_I_Ca_Na = []
        self.list_I_Ca_K  = []
        self.list_I_Ca_t  = []

        self.list_d       = []
        self.list_f       = []
        #3.I_K
        self.list_alpha_X = []
        self.list_beta_X  = []
        self.list_X       = []
        self.list_Xi      = []
        self.list_I_K     = []
        #4I_K1
        self.list_alpha_K1= []
        self.list_beta_K1 = []
        self.list_K1      = []
        self.list_I_K1    = []
        self.list_I_K1_inf  = []
        #5I_Kp
        self.list_Kp     = []
        self.list_I_Kp    = []
        #6I_EX_Na_Ca
        self.list_I_EX_Na_Ca   = []
        #6list_Ib
        self.list_Ib     = []
        #I_ion
        self.list_I_ion   = []
        #I_test
        self.list_alpha_test = []
        self.list_beta_test  = []
        self.list_test       = []
        self.list_I_test    = []

    def plot(self,para_a,para_b,length):
        #this is to randomly plot a figure.
        import numpy as np
        import matplotlib.pyplot as plt
        p = plt.subplot(111)
        p.plot(para_a[0:length],para_b[0:length],color = 'black',linewidth = 0.5)
        p.spines['right'].set_visible(False)
        p.spines['top'].set_visible(False)
        p.yaxis.set_ticks_position('left')
        p.xaxis.set_ticks_position('bottom')
        plt.ylabel('Voltage (mV)')
        plt.xlabel('Time (s)')
        plt.axis([5 , 100, -100, 60])
        plt.xticks(fontsize = 32)
        plt.yticks(fontsize = 32)
        #plt.savefig(name[0:7] + ' total.pdf')
        plt.show()        

    def init_V(self,V):
        self.list_V.append(V)

    def init_Na(self):    
        #____(alpha, beta) of (m,h,j)____________________________________
        self.list_alpha_m.append(alpha_m(self.list_V[-1]))
        self.list_beta_m.append(beta_m(self.list_V[-1]))
        self.list_alpha_h.append(alpha_h(self.list_V[-1])) 
        self.list_beta_h.append(beta_h(self.list_V[-1]))
        self.list_alpha_j.append(alpha_j(self.list_V[-1]))
        self.list_beta_j.append(beta_j(self.list_V[-1]))
        #____(m,h,j) inf as initiation for calculating the current_______
        self.list_m.append(probability_inf(self.list_alpha_m[-1],self.list_beta_m[-1]))
        self.list_h.append(probability_inf(self.list_alpha_h[-1],self.list_beta_h[-1]))
        self.list_j.append(probability_inf(self.list_alpha_j[-1],self.list_beta_j[-1]))
        self.list_I_Na.append(gbar_Na *  self.list_m[-1]**3 * self.list_h[-1] * self.list_j[-1] * (self.list_V[-1]-E_Na))
        #########################____________________finish
        
    def init_Ca(self):
        #I_Ca_for ________d,f_______________________________
        '''
        previously, the steps are:
            1. calculate the alpha, beta of GATE
            2. calculate the INF_GATE as the GATE[0]
            3. accumulation GATEs based on the integral.

        Now it should be:
            1. calculate the INF_GATE as the 1st GATE
            2. calculate the alpha, beta of the GATE base on the INF_GATE
            3. accumulation GATEs based on the integral.
        '''
        #____(d,f,f_Ca)_______________________________
        self.list_inf_d.append(inf_d(self.list_V[-1]))
        self.list_inf_f.append(inf_f(self.list_V[-1]))

        self.list_tau_d.append(tau_d(self.list_V[-1]))
        self.list_tau_f.append(tau_f(self.list_V[-1]))
        
        self.list_alpha_d.append(self.list_inf_d[-1]/self.list_tau_d[-1])
        self.list_beta_d.append((1-self.list_inf_d[-1])/self.list_tau_d[-1])
        self.list_alpha_f.append(self.list_inf_f[-1]/self.list_tau_f[-1])
        self.list_beta_f.append((1-self.list_inf_f[-1])/self.list_tau_f[-1])
        
        self.list_d.append(inf_d(self.list_V[-1]))
        self.list_f.append(inf_f(self.list_V[-1]))
        
        self.list_f_Ca.append(f_Ca(self.list_V[-1]))

        I_Ca_Ca = I_general((self.list_V[-1]), P_Ca_Ca, 2 ,Gama_Cai, Gama_Cao, Ca_i, Ca_o)
        self.list_I_Ca_Ca.append(I_Ca_Ca) 
        I_Ca_Na = I_general((self.list_V[-1]), P_Ca_Na, 1 ,Gama_Nai, Gama_Nao, Na_i, Na_o)
        self.list_I_Ca_Na.append(I_Ca_Na) 
        I_Ca_K  = I_general((self.list_V[-1]), P_Ca_K , 1 ,Gama_Ki , Gama_Ko , K_i , K_o )
        self.list_I_Ca_K.append(I_Ca_K)
        I_Ca_t = I_Ca_Ca + I_Ca_Na + I_Ca_K
        self.list_I_Ca_t.append(  self.list_d[-1]    * self.list_f[-1]  * self.list_f_Ca[-1]    *    I_Ca_t        )     
        #########################____________________finish
        
    def init_K(self):  
        #I_K_for potassium________X_______________________________
        self.list_alpha_X.append(alpha_X(self.list_V[-1]))
        self.list_beta_X.append(beta_X(self.list_V[-1]))
        
        self.list_X.append(self.list_alpha_X[-1] /(self.list_alpha_X[-1]+self.list_beta_X[-1]))
        self.list_Xi.append(Xi(self.list_V[-1]))
        self.list_I_K.append (gbar_K  *  self.list_X[-1]    * self.list_Xi[-1]**2 *                   (self.list_V[-1]-E_K ))        
        #########################____________________finish
    def init_K1(self):  
        #I_K1 for potassium_______K1__________________________
        self.list_alpha_K1.append(alpha_K1(self.list_V[-1]))
        self.list_beta_K1.append(beta_K1(self.list_V[-1]))
        self.list_K1.append(self.list_alpha_K1[-1] /(self.list_alpha_K1[-1]+self.list_beta_K1[-1]))
        self.list_I_K1.append(gbar_K1  *  self.list_K1[-1]               *                   (self.list_V[-1]-E_K1))
        #########################____________________finish
        
    def init_Kp(self):  
        #I_Kp for potassium
        self.list_Kp.append(Kp(self.list_V[-1]))
        self.list_I_Kp.append(gbar_Kp  * self.list_Kp[-1]   *                                (self.list_V[-1]-E_Kp))

    def init_EX_Na_Ca(self):
        #I_EX_Na_Ca
        self.list_I_EX_Na_Ca.append(I_EX_Na_Ca(self.list_V[-1]))        
    def init_Ib(self):
        self.list_Ib.append(Ib(self.list_V[-1]))        

    def init_I_ion(self):
        self.list_I_ion.append(I[i] -self.list_I_EX_Na_Ca[-1] - self.list_I_Na[-1]- self.list_I_Ca_t[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_Ib[-1])        

    def init(self,V):
        self.init_V(V)
        self.init_Na()
        self.init_Ca()
        self.init_K()
        self.init_K1()
        self.init_Kp()
        self.init_EX_Na_Ca()
        self.init_Ib()
        self.init_I_ion()

#then we got all the values for the first round, then here comes the integral
#since the delta_t is 0.01, so the first thing that change is V
#then the V cause another changes.
        
    def alter_V(self):
        self.list_V.append(self.list_V[-1] +  deltaT * self.list_I_ion[-1])
        
    def alter_Na(self):        
        self.list_alpha_m.append(alpha_m(self.list_V[-1]))
        self.list_beta_m.append(beta_m(self.list_V[-1]))
        self.list_alpha_h.append(alpha_h(self.list_V[-1]))
        self.list_beta_h.append(beta_h(self.list_V[-1]))
        self.list_alpha_j.append(alpha_j(self.list_V[-1]))
        self.list_beta_j.append(beta_j(self.list_V[-1]))
        
        self.list_m.append( self.list_m[-1]+ deltaT * (self.list_alpha_m[-1]* (1-self.list_m[-1])- self.list_beta_m[-1]* self.list_m[-1]))
        self.list_h.append(self.list_h[-1]+ deltaT * (self.list_alpha_h[-1]* (1-self.list_h[-1])- self.list_beta_h[-1]* self.list_h[-1]))
        self.list_j.append(self.list_j[-1]+ deltaT * (self.list_alpha_j[-1]* (1-self.list_j[-1])- self.list_beta_j[-1]* self.list_j[-1]))
        
        self.list_I_Na.append(gbar_Na *  self.list_m[-1]**3 * self.list_h[-1] * self.list_j[-1]*(self.list_V[-1]-E_Na))         


    def alter_Ca_t(self):

        
        self.list_inf_d.append(inf_d(self.list_V[-1]))
        self.list_inf_f.append(inf_f(self.list_V[-1]))
        self.list_tau_d.append(tau_d(self.list_V[-1]))
        self.list_tau_f.append(tau_f(self.list_V[-1]))

        self.list_alpha_d.append(self.list_inf_d[-1]/self.list_tau_d[-1])
        self.list_beta_d.append((1-self.list_inf_d[-1])/self.list_tau_d[-1])
        self.list_alpha_f.append(self.list_inf_f[-1]/self.list_tau_f[-1])
        self.list_beta_f.append((1-self.list_inf_f[-1])/self.list_tau_f[-1])


        self.list_f_Ca.append(f_Ca(self.list_V[-1]))


        self.list_d.append( self.list_d[-1]+ deltaT * (self.list_alpha_d[-1]* (1-self.list_d[-1])- self.list_beta_d[-1]* self.list_d[-1]))
        self.list_f.append(self.list_f[-1]+ deltaT * (self.list_alpha_f[-1]* (1-self.list_f[-1])- self.list_beta_f[-1]* self.list_f[-1]))
       


        I_Ca_Ca = I_general((self.list_V[-1]),P_Ca_Ca,2,Gama_Cai,Gama_Cao,Ca_i,Ca_o)
        self.list_I_Ca_Ca.append(I_Ca_Ca) 
        I_Ca_Na = I_general((self.list_V[-1]),P_Ca_Na,1,Gama_Nai,Gama_Nao,Na_i,Na_o)
        self.list_I_Ca_Na.append(I_Ca_Na) 
        I_Ca_K  = I_general((self.list_V[-1]),P_Ca_K ,1,Gama_Ki,Gama_Ko,K_i,K_o)
        self.list_I_Ca_K.append(I_Ca_K) 

        I_Ca_t = I_Ca_Ca + I_Ca_Na + I_Ca_K

        self.list_I_Ca_t.append(  self.list_d[-1]    * self.list_f[-1]  * self.list_f_Ca[-1]    *    I_Ca_t        )
        


    def alter_I_K(self):
        self.list_alpha_X.append(alpha_X(self.list_V[-1]))
        self.list_beta_X.append(beta_X(self.list_V[-1]))

        self.list_X.append(self.list_X[-1]+ deltaT * (self.list_alpha_X[-1]* (1-self.list_X[-1])- self.list_beta_X[-1]* self.list_X[-1]))
        self.list_Xi.append(Xi(self.list_V[-1]))
        
        self.list_I_K.append (gbar_K  *  self.list_X[-1]    * self.list_Xi[-1]**2 *                   (self.list_V[-1]-E_K ))

    def alter_I_K1(self):
        self.list_alpha_K1.append(alpha_K1(self.list_V[-1]))
        self.list_beta_K1.append(beta_K1(self.list_V[-1]))
        self.list_K1.append(self.list_alpha_K1[-1] /(self.list_alpha_K1[-1]+self.list_beta_K1[-1]))
        self.list_I_K1.append(gbar_K1  *  self.list_K1[-1]               *                   (self.list_V[-1]-E_K1))
    def alter_I_Kp(self):
        self.list_Kp.append(Kp(self.list_V[-1]))
        self.list_I_Kp.append(gbar_Kp  * self.list_Kp[-1]   *                                (self.list_V[-1]-E_Kp))
    def alter_Ib(self):
        self.list_Ib.append(Ib(self.list_V[-1]))
    def alter(self,n):
        for i in range(n):
            self.alter_V() 
            self.alter_Na()
            self.alter_I_si()
            self.alter_I_K()
            self.alter_I_K1()
            self.alter_I_Kp()
            self.alter_Ib()
            self.list_I_ion.append(I[i] - self.list_I_Na[-1]- self.list_I_si[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_Ib[-1])

        
    def optionalter(self,n,a = 1,b = 1,c = 1,d = 1, e = 1, f = 1):
        for i in range(n):
            self.alter_V()
            if a == 1:
                self.alter_Na()
            else:
                self.list_I_Na.append(0)
                
            if b == 1:
                self.alter_Ca_t()
            else:
                self.list_I_Ca_t.append(0)
                
            if c == 1:
                self.alter_I_K()
            else:
                self.list_I_K.append(0)
                
            if d == 1:
                self.alter_I_K1()
            else:
                self.list_I_K1.append(0)
                
            if e == 1:
                self.alter_I_Kp()
            else:
                self.list_I_Kp.append(0)
                
            if f == 1:
                self.alter_Ib()
            else:
                self.list_Ib.append(0)
  
            self.list_I_ion.append(I[i] - self.list_I_Na[-1]- self.list_I_Ca_t[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_Ib[-1])        


a = actionpo()

zzz = 500000
a.init(-85) 



a.optionalter(zzz,1,1,1,1,1,0)
'''
try:

    
    a.optionalter(zzz,1,0,1,1,1,1)
except:
    alllist = [a.list_V,
           a.list_alpha_m, a.list_beta_m, a.list_alpha_h, a.list_beta_h, a.list_alpha_j, a.list_beta_j,
           a.list_h, a.list_h, a.list_j,a.list_I_Na,
           a.list_alpha_d, a.list_beta_d, a.list_alpha_f, a.list_beta_f,
           a.list_d, a.list_f,a.list_I_si,
           a.list_alpha_X, a.list_beta_X, 
           a.list_X, a.list_Xi, a.list_I_K,
           a.list_alpha_K1, a.list_beta_K1,  
           a.list_K1, a.list_I_K1,  
           a.list_Kp, a.list_I_Kp,
           a.list_Ib,
           a.list_I_ion,
           ]
'''








p = plt.subplot(111)



p.plot(t[0:len(a.list_V)],a.list_V[0:len(a.list_V)],color = 'k',linewidth = 3)
p.plot(t[0:len(a.list_V)],a.list_I_Ca_t[0:len(a.list_V)],color = 'r',linewidth = 3)
#p.plot(t[0:zzz],a.list_I_Na[0:zzz],color = 'r',linewidth = 1.5)
 
#p.plot(t[0:zzz],a.list_V[0:zzz],color = 'g',linewidth = 1.5)
#p.plot(t[0:zzz],a.list_h[0:zzz],color = 'k',linewidth = 1.5)
#p.plot(t[0:zzz],a.list_I_Na[0:zzz],color = 'red',linewidth = 1.5)
#p.plot(t[0:50000],a.list_I_K [0:50000],color = 'blue',linewidth = 1.5)
#p.plot(t[0:10000],list_I_L[0:10000],color = 'green',linewidth = 1.5)
#p.plot(a.list_V[0:zzz],a.list_I_Kp[0:zzz],color = 'g',linewidth = 1.5)

p.spines['right'].set_visible(False)
p.spines['top'].set_visible(False)
p.yaxis.set_ticks_position('left')
p.xaxis.set_ticks_position('bottom')

p.text(320,-52,'stim', fontsize = 28,color = 'green')
p.text(1000,-45,"4 folds'", fontsize = 28, color = 'green')
p.text(1050,-55,"stim", fontsize = 28, color = 'green')
p.arrow(450, -50, 30, 0, head_width=5, head_length=5, fc='g', ec='g')
p.arrow(1200, -50, 30, 0, head_width=5, head_length=5, fc='g', ec='g')
p.arrow(1400, 15, 20, -10, head_width=5, head_length=5, fc='r', ec='r')
p.arrow(700, 20, -50, -10, head_width=5, head_length=5, fc='k', ec='k')
p.text(1400,20,"calcium conductance", fontsize = 28, color = 'r')
p.text(700,20,"action potentials", fontsize = 28)
plt.xticks(fontsize = 32)
plt.yticks(fontsize = 32)

plt.ylabel('Voltage (mV)',fontsize = 28)
plt.xlabel('Time (mS)', fontsize = 28)
plt.axis([200 ,2000, -100, 80])
plt.show()








