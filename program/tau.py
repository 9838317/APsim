#this is for import all the math fuctions from another file, and confirm that the program could be working under different OS.
import os
import sys
current_path = os.getcwd()                              #get current dir
parent_path = os.path.dirname(current_path)             #get parent  dir
sys.path.append(current_path)                           
sys.path.append(parent_path)                            #to append the previous

#import constant
from z_pa_and_con import *
import numpy as np
import matplotlib.pyplot as plt

def probability_inf(a,b):
    return 1.0 * a / (a + b)

def tau(a,b):
    return 1.0 / (a + b)

'''
This is a program use the formular of:
m_0 = m_inf = alpha_m / (alpha_m + beta_m)
and dy/dt = (y_inf - y) / tau
'''

#first cycle:
#CALCULATE alpha,beta, then m,n
class actionpo(object):
    def __init__(self):
        #Create all the list
        #0.V        
        self.list_V              = []
        #1.I_Na
        self.list_alpha_m        = []
        self.list_beta_m         = []
        self.list_alpha_h        = []
        self.list_beta_h         = []
        self.list_alpha_j        = []
        self.list_beta_j         = []
        self.list_m              = []
        self.list_h              = []
        self.list_j              = []
        self.list_I_Na           = []
        #2.I_si
        self.list_alpha_d        = []
        self.list_beta_d         = []
        self.list_alpha_f        = []
        self.list_beta_f         = []
        self.list_d              = []
        self.list_f              = []
        self.list_I_si           = []
        #3.I_K
        self.list_alpha_X        = []
        self.list_beta_X         = []
        self.list_X              = []
        self.list_Xi             = []
        self.list_I_K            = []
        #4.I_K1
        self.list_alpha_K1       = []
        self.list_beta_K1        = []
        self.list_K1             = []
        self.list_I_K1           = []
        self.list_I_K1_inf       = []
        #5.I_Kp
        self.list_Kp             = []
        self.list_I_Kp           = []
        #6.I_b
        self.list_I_b            = []
        #I_ion
        self.list_I_ion          = []
        #I_test
        self.list_alpha_test     = []
        self.list_beta_test      = []
        self.list_test           = []
        self.list_I_test         = []
        #########################____________________finish
        
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
        #plt.axis([-1 , 4 0, -100, 60])
        #plt.savefig(name[0:7] + ' total.pdf')
        plt.show()        

    def init_V(self,V):
        self.list_V.append(V)
        #########################____________________finish
        
    def init_I_Na(self):    
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
        
    def init_I_si(self):  
        #____(alpha, beta) of (d,f)_______________________________
        self.list_alpha_d.append(alpha_d(self.list_V[-1]))
        self.list_beta_d.append(beta_d(self.list_V[-1]))
        self.list_alpha_f.append(alpha_f(self.list_V[-1]))
        self.list_beta_f.append(beta_f(self.list_V[-1]))
        #____(d,f) inf as initiation for calculating the current_______
        self.list_d.append(probability_inf(self.list_alpha_d[-1],self.list_beta_d[-1]))
        self.list_f.append(probability_inf(self.list_alpha_f[-1],self.list_beta_f[-1]))
        self.list_I_si.append(gbar_si *  self.list_d[-1]    * self.list_f[-1]  *                (self.list_V[-1]-E_si))
        #########################____________________finish
        
    def init_I_K(self):  
        #____(alpha, beta) of (X)_______________________________
        self.list_alpha_X.append(alpha_X(self.list_V[-1]))
        self.list_beta_X.append(beta_X(self.list_V[-1]))
        #____(X) inf as initiation for calculating the current_______
        self.list_X.append(probability_inf(self.list_alpha_X[-1],self.list_beta_X[-1]))
        self.list_Xi.append(Xi(self.list_V[-1]))
        self.list_I_K.append (gbar_K  *  self.list_X[-1]    * self.list_Xi[-1] *                   (self.list_V[-1]-E_K ))
        #########################____________________finish
        
    def init_I_K1(self):  
        #I_K1 for potassium_______K1__________________________
        self.list_alpha_K1.append(alpha_K1(self.list_V[-1]))
        self.list_beta_K1.append(beta_K1(self.list_V[-1]))
        self.list_K1.append(probability_inf(self.list_alpha_K1[-1],self.list_beta_K1[-1]))
        self.list_I_K1.append(gbar_K1  *  self.list_K1[-1]               *                   (self.list_V[-1]-E_K1))
        #########################____________________finish
        
    def init_I_Kp(self):  
        #I_Kp for potassium
        self.list_Kp.append(Kp(self.list_V[-1]))
        self.list_I_Kp.append(gbar_Kp  * self.list_Kp[-1]   *                                (self.list_V[-1]-E_Kp)) 
        #########################____________________finish
        
    def init_I_b(self):
        self.list_I_b.append(I_b(self.list_V[-1]))
        #########################____________________finish

    def init_I_ion(self):  
        self.list_I_ion.append(I[i] - self.list_I_Na[-1]- self.list_I_si[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_I_b[-1])        
        #########################____________________finish

    def init(self,V):
        self.init_V(V)
        self.init_I_Na()
        self.init_I_si()
        self.init_I_K()
        self.init_I_K1()
        self.init_I_Kp()
        self.init_I_b()
        self.init_I_ion()
        #########################____________________finish
        
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
        temp_1 = probability_inf(self.list_alpha_m[-1],self.list_beta_m[-1])
        temp_2 = probability_inf(self.list_alpha_h[-1],self.list_beta_h[-1])

        temp_3 = probability_inf(self.list_alpha_j[-1],self.list_beta_j[-1])

        self.list_m.append(self.list_m[-1]+ deltaT * ((temp_1 - self.list_m[-1])* (self.list_alpha_m[-1]+self.list_beta_m[-1])))
        self.list_h.append(self.list_h[-1]+ deltaT * ((temp_2 - self.list_h[-1])* (self.list_alpha_h[-1]+self.list_beta_h[-1])))
        self.list_j.append(self.list_j[-1]+ deltaT * ((temp_3 - self.list_j[-1])* (self.list_alpha_j[-1]+self.list_beta_j[-1])))
        self.list_I_Na.append(gbar_Na *  self.list_m[-1]**3 * self.list_h[-1] * self.list_j[-1]*(self.list_V[-1]-E_Na))
  

        
    def alter_I_si(self):        
        self.list_alpha_d.append(alpha_d(self.list_V[-1]))
        self.list_beta_d.append(beta_d(self.list_V[-1]))
        self.list_alpha_f.append(alpha_f(self.list_V[-1]))
        self.list_beta_f.append(beta_f(self.list_V[-1]))
        temp_1 = probability_inf(self.list_alpha_d[-1],self.list_beta_d[-1])
        temp_2 = probability_inf(self.list_alpha_f[-1],self.list_beta_f[-1])
        self.list_d.append(self.list_d[-1]+ deltaT * ((temp_1 - self.list_d[-1])* (self.list_alpha_d[-1]+self.list_beta_d[-1])))
        self.list_f.append(self.list_f[-1]+ deltaT * ((temp_2 - self.list_f[-1])* (self.list_alpha_f[-1]+self.list_beta_f[-1])))
        self.list_I_si.append(gbar_si *  self.list_d[-1]    * self.list_f[-1]  *                (self.list_V[-1]-E_si))
    def alter_I_K(self):
        self.list_alpha_X.append(alpha_X(self.list_V[-1]))
        self.list_beta_X.append(beta_X(self.list_V[-1]))
        temp_1 = probability_inf(self.list_alpha_X[-1],self.list_beta_X[-1])
        self.list_X.append(self.list_X[-1]+ deltaT * ((temp_1 - self.list_X[-1])* (self.list_alpha_X[-1]+self.list_beta_X[-1])))
        self.list_Xi.append(Xi(self.list_V[-1]))
        self.list_I_K.append (gbar_K  *  self.list_X[-1]    * self.list_Xi[-1] *                   (self.list_V[-1]-E_K ))

    def alter_I_K1(self):
        self.list_alpha_K1.append(alpha_K1(self.list_V[-1]))
        self.list_beta_K1.append(beta_K1(self.list_V[-1]))
        self.list_K1.append(probability_inf(self.list_alpha_K1[-1],self.list_beta_K1[-1]))
        self.list_I_K1.append(gbar_K1  *  self.list_K1[-1]               *                   (self.list_V[-1]-E_K1))
    def alter_I_Kp(self):
        self.list_Kp.append(Kp(self.list_V[-1]))
        self.list_I_Kp.append(gbar_Kp  * self.list_Kp[-1]   *                                (self.list_V[-1]-E_Kp))
    def alter_I_b(self):
        self.list_I_b.append(I_b(self.list_V[-1]))
    def alter(self,n):
        for i in range(n):
            self.alter_V() 
            self.alter_Na()
            self.alter_I_si()
            self.alter_I_K()
            self.alter_I_K1()
            self.alter_I_Kp()
            self.alter_I_b()
            self.list_I_ion.append(I[i] - self.list_I_Na[-1]- self.list_I_si[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_I_b[-1])
    def optionalter(self,n,a = 1,b = 1,c = 1,d = 1, e = 1, f = 1):
        for i in range(n):
            self.alter_V()
            if a == 1:
                self.alter_Na()
            else:
                self.list_I_Na.append(0)
                
            if b == 1:
                self.alter_I_si()
            else:
                self.list_I_si.append(0)
                
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
                self.alter_I_b()
            else:
                self.list_I_b.append(0)
  
            self.list_I_ion.append(I[i] - self.list_I_Na[-1]- self.list_I_si[-1]-self.list_I_K[-1] - self.list_I_K1[-1] -  self.list_I_Kp[-1] - self.list_I_b[-1])        

    def test(self, k):
        self.list_I_test = [self.list_I_ion, self.list_I_Na,self.list_I_si,self.list_I_K, self.list_I_K1, self.list_I_Kp, self.list_I_b ]
        for i in self.list_I_test:
            print i[k]
        print '----'
        self.list_test = [self.list_alpha_m,self.list_beta_m,self.list_alpha_h,self.list_beta_h,self.list_alpha_j,self.list_beta_j ]
        for i in self.list_test:
            print i[k]  
        print '----'
        self.list_integral = [self.list_m,self.list_h,self.list_j]
        for i in self.list_integral:
            print i[k]  
        print '----'      

a = actionpo()

datapoint = 100000
a.init(-85) 


a.optionalter(datapoint,1,1,1,1,1,1)
a.plot(t,a.list_V,datapoint)
a.test(1)
'''
try:
    a.optionalter(zzz,1,1,1,1,1,1)
except:
    alllist = [a.list_V,
           a.list_alpha_m, a.list_beta_m, a.list_alpha_h, a.list_beta_h, a.list_alpha_j, a.list_beta_j,
           a.list_h, a.list_h, a.list_j,a.list_I_Na,
           a.list_alpha_d, a.list_beta_ d, a.list_alpha_f, a.list_beta_f,
           a.list_d, a.list_f,a.list_I_si,
           a.list_alpha_X, a.list_beta_X, 
           a.list_X, a.list_Xi, a.list_I_K,
           a.list_alpha_K1, a.list_beta_K1,  
           a.list_K1, a.list_I_K1,  
           a.list_Kp, a.list_I_Kp,
           a.list_I_b,
           a.list_I_ion,
           ]
'''






'''
p = plt.subplot(111)


 
p.plot(t[0:len(a.list_V)],a.list_V[0:len(a.list_V)],color = 'k',linewidth = 1.5)

p = plt.subplot(111)
p.plot(t[0:zzz],a.list_I_Na[0:zzz],color = 'r',linewidth = 1.5)
 
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
#p.text(30,-35,'Bigger')
#p.text(35,-45,'stim')
#p.text(200,45,'cell with high calcium current')
#p.text(200,35,'sodium current = 0%')
p.arrow(45, -50, 30, 0, head_width=5, head_length=5, fc='k', ec='k')


plt.ylabel('Voltage (mV)')
plt.xlabel('Time (mS)')
#plt.axis([0 ,800, -100, 80])
plt.show()

'''






