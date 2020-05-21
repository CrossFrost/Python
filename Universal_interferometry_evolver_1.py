# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:58:55 2020
Only solves for discrete k
@author: Jay
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import special 

#initializing g,t 

h = 1e-9 #step size in t in seconds 
t = np.arange(0,1e-6+h,h) #time array
grid_size_g = 31 #always put odd integer here
size_g = (grid_size_g,len(t)) #size of the evolution array
size_rk_factors = (grid_size_g,4) #size of RK coefficient array
K = np.zeros(size_rk_factors,dtype=complex) #array for storing RK coefficients
g = np.zeros(size_g,dtype=complex) #evolution data stored here

#set initial condition
g[np.int((grid_size_g-1)*0.5),0] = 1 

#define derivatives 
def dxdt(t,g):
    
    dxdt = np.zeros(grid_size_g,dtype=complex) #initiate coefficients 
    recoil = 2*np.pi*3.8*1e3 #for Rb87 and 780 nm light 
    omega  = 80*recoil #Rabi frequency 
    for i in range(len(g)-4): #first two and last two elements of the array won't evolve
        dxdt[i+2] = -1j*((recoil*(i+2-(grid_size_g-1)*0.5)**2 + omega)*g[i+2] + 0.5*omega*(g[i+4]+g[i]))
    return dxdt    

for i in range(len(t)-1): #loop RK4 
    
    temp = np.zeros(grid_size_g,dtype=complex) #temporary variable 

    K[:,0] = dxdt(t[i],g[:,i]) #Runge Kutta coefficients
    
    temp = g[:,i]+0.5*K[:,0]*h

    K[:,1] = dxdt(t[i]+h/2.0,temp)
    
    temp = g[:,i]+0.5*K[:,1]*h

    K[:,2] = dxdt(t[i]+h/2.0,temp)

    temp = g[:,i]+K[:,2]*h

    K[:,3] = dxdt(t[i]+h,temp)

    g[:,i+1] = g[:,i] + h*(K[:,0]+2*K[:,1]+2*K[:,2]+K[:,3])/6.0 #evolver step 

#create Bessel grid 
N = np.int((grid_size_g-1)*0.5)
Bessel_array = np.zeros(grid_size_g)
phi_d = 1.9 #Kick strength 1.9 for omega = 80*recoil 
for i in range(grid_size_g):
    if (i+1)%2 == 0:
        Bessel_array[i] = special.jv(0.5*(i-N),phi_d)**2
    else : 
        Bessel_array[i] = 0
            
plt.plot(np.arange(-N,N+1,1),np.abs(g[:,len(t)-1])**2,label='Simulations') #plots distribution at the end of evolution 
plt.plot(np.arange(-N,N+1,1),Bessel_array,label='Bessel function',linestyle='--',color='r')
plt.xlabel("Momentum orders (n h_bar k)")
plt.ylabel("Probability")
plt.legend(loc='upper right')
plt.show






    
    