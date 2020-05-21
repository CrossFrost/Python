# -*- coding: utf-8 -*-
"""
Created on Thu May 21 10:24:14 2020

@author: Jay
"""


import matplotlib.pyplot as plt
import numpy as np

#initializing g,t

h               = 50e-9 #step size in t in seconds 
t_end           = 3e-6 #total evolution period in seconds
t               = np.arange(0,t_end+h,h) #time array 
k               = 30 #number of discrete orders -k to k 
beta            = 101 #quasimomentum distribution within an order -0.5 to 0.5 
grid_size_g     = (2*k+1)*beta 
size_g          = (grid_size_g,len(t)) #size of the evolution array
size_rk_factors = (grid_size_g,4) #size of RK coefficient array
K               = np.zeros(size_rk_factors,dtype=complex) #array for storing RK coefficients
g               = np.zeros(size_g,dtype=complex) #evolution coefficient stored here

#set initial condition 

def func(x,sigma,mu): #Gaussian distribution function for initial distribution
    return (1/(sigma*np.sqrt(2*np.pi)))*(np.exp(-0.5*((x-mu)/sigma)**2))

p      = np.arange(0,grid_size_g,1) 
sigma  = 0.25 #in units of 2 h_bar k, k=2*pi/lambda  
mu     = np.float(50+k*beta) #zero of the momentum array
g[:,0] = func(p,sigma*2*beta,mu)
Norm   = (np.abs(np.sum(g[:,0]**2)))**0.5 #normalization constant so that sum(g**2)=1 
g[:,0] = (1.0/Norm)*g[:,0]  #Gaussian initial distribution realized 

#g[50+k*beta,0] = 1 

#define derivatives 

def dxdt(t,g):
    
    dxdt   = np.zeros(grid_size_g,dtype=complex) #initiate coefficients 
    recoil = 2*np.pi*3.8*1e3 #for Rb87 and 780 nm light 
    omega  = 50*recoil #Rabi frequency 
    x      = np.zeros(len(g))
    
    for i in range(len(x)):   #set up the index array 
        delta = (np.mod(i,beta)-(beta-1)*0.5)*(1.0/(beta-1)) #quasimomentum index 
        I     = i/beta #discrete momentum index
        x[i]  = delta + I - k #index array 
        
    for i in range(len(g)-4*beta): #first and last two I-indicies do not evolve 
        dxdt[i+2*beta] = -1j*((recoil*(x[i+2*beta])**2 + omega)*g[i+2*beta] + 0.5*omega*(g[i+4*beta] + g[i])) 
    return dxdt     
        
for i in range(len(t)-1): #loop RK4 
    
    temp     = np.zeros(grid_size_g,dtype=complex) #temporary variable 
    K[:,0]   = dxdt(t[i],g[:,i]) #Runge Kutta coefficients     
    temp     = g[:,i]+0.5*K[:,0]*h
    K[:,1]   = dxdt(t[i]+h/2.0,temp)   
    temp     = g[:,i]+0.5*K[:,1]*h
    K[:,2]   = dxdt(t[i]+h/2.0,temp)
    temp     = g[:,i]+K[:,2]*h
    K[:,3]   = dxdt(t[i]+h,temp)
    g[:,i+1] = g[:,i] + h*(K[:,0]+2*K[:,1]+2*K[:,2]+K[:,3])/6.0 #evolver step         
    
    if i%(0.1*(len(t)-1)) == 0: #displays % completion of loop in units of 10%
        q = i/(0.1*(len(t)-1)) 
        print str((q+1)*10) + '% complete'

#setup x-axis for plotting
x_plot = np.zeros(grid_size_g) 
for i in range(len(x_plot)):  
    delta = (np.mod(i,beta)-(beta-1)*0.5)*(1.0/(beta-1))
    I     = i/beta
    x_plot[i]  = delta + I - k

#plots momentum distribution at the end of evolution 
plt.plot(x_plot,np.abs(g[:,len(t)-1])**2)    
plt.xlabel("Momentum orders (n h_bar k)")
plt.ylabel("Probability")    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
