# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:32:54 2020

@author: Jay
"""


import numpy as np
from scipy import special 
import matplotlib.pyplot as plt

data = np.loadtxt("data_for_phase_asymm.txt",unpack =True) #experimental data import

elle = 2 #scaled time between the kicks, elle = 2 for resonant time (default = 2)
sigma = 0.03 #scaled width of initial distribution (default = 0.03)
mu = 0 #offset from zero in momentum space
phi_d = 1.15 #kick strength (default = 1.15)
#phase = np.pi - np.pi*0.3

def func(x,sigma,mu): #Gaussian distribution function for initial distribution
    return (1/(sigma*np.sqrt(2*np.pi)))*(np.exp(-0.5*((x-mu)/sigma)**2))

def func_2(x,phi_d,elle,phase,order): #Bessel distribution function for distributuon modified after kicks
    return (special.jv(order,(2*phi_d*np.cos((np.pi*(1+2*x)*elle-phase)/2.0))))**2 

x_array = np.arange(-0.5,0.5,0.001) #basis array of quasimomenta

def asymm_func(phase,order): #gives measure of asymmetry as an output for a given lattice phase shift

    pd_array = np.zeros(np.size(x_array)-1) #stores resultant distribution
    
    for i in range(np.size(pd_array)): #core formula
        pd_array[i] = func_2(x_array[i],phi_d,elle,phase,order)*func(x_array[i],sigma,mu)*0.001
        
    #l1 = np.int(np.floor(np.size(pd_array)/2.0)) #slicing distribution to calculate asymmetry
    #l2 = np.int(np.ceil(np.size(pd_array)/2.0)) 
        
    #asymmetry = -(np.sum(pd_array[0:l1])-np.sum(pd_array[l2:np.size(pd_array)-1]))/(np.sum(pd_array))
    return pd_array 

#phase_array = np.arange(0,2*np.pi+2*np.pi/100,2*np.pi/100) #important to include 2*pi as last entry

number_of_orders = 1
size = (2*number_of_orders+1,(np.size(x_array)-1))
distribution = np.zeros(size)

#for i in range(np.size(phase_array)): #fills up asymmetry array 

s_temp = np.zeros(2)

for i in range(3):    
    distribution[i] = asymm_func(np.pi,i-1) 
    #plt.plot(x_array[0:np.size(x_array)-1],distribution[i])
    s_temp = np.concatenate((s_temp,distribution[i]))

s_temp = np.delete(s_temp,[0,1])    

s = np.zeros(2)

for i in range(3):    
    distribution[i] = asymm_func(60*(np.pi/180),i-1) 
    #plt.plot(x_array[0:np.size(x_array)-1],distribution[i])
    s = np.concatenate((s,distribution[i]))

s = np.delete(s,[0,1]) 

s_1 = np.zeros(2)

for i in range(3):    
    distribution[i] = asymm_func(250*(np.pi/180),i-1) #280 seems the right place
    #plt.plot(x_array[0:np.size(x_array)-1],distribution[i])
    s_1 = np.concatenate((s_1,distribution[i]))

s_1 = np.delete(s_1,[0,1]) 
 
x_array_2 = np.arange(-1.5,1.5,3.0/(1*np.size(s)))
plt.axvline(x=-1,linestyle='--',color='r')
plt.axvline(x=-0,linestyle='--',color='r')
plt.axvline(x=1,linestyle='--',color='r')
#plt.plot(x_array_2[0+200:np.size(x_array_2)-200],s_temp[0+200:np.size(x_array_2)-200])
#plt.plot(x_array_2,s_temp)
plt.plot(x_array_2[0:np.size(x_array_2)],s_1[0:np.size(x_array_2)])
plt.show()
#plt.plot(phase_array,asymmetry_array)

#plt.errorbar(data[0]*(np.pi/180), data[1],linestyle='none',marker='o',color='r', yerr=data[2],label='experiment')
#plt.plot(phase_array,asymmetry_array,label='theory')
#plt.xlabel("Phase offset (radians)")
#plt.ylabel("Normalized asymmetry")
#plt.legend(loc='upper left')
#plt.savefig('asymmetry_data_plot_1.svg')
#plt.show()
#print asymmetry