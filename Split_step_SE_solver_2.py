# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:14:17 2020

@author: Jay
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy.fftpack import fft,ifft

#INPUT PARAMETERS

#set time parameters
t0 = 0
tf = 40000
dt = 0.1

#set space grid parameters
dx = 0.1
xmax = 8
#xmin will be -xmax, making the situation symmetric

#how far seperated are the potential minima
seperation = 6
#minima at seperation/2 and -seperation/2

#How often we measure occupation and measure state .
nmeasure = 100

t = np.arange(t0,tf+dt,dt)
x = np.arange(-xmax,xmax+dx,dx)

nt = len(t)
N = len(x)

k = np.concatenate((np.arange(0,(N-1)*0.5+1,1),np.arange(-80,0,1)),axis=0)*(2*np.pi/(N*dx))

size = (np.int(np.floor(nt/nmeasure)),2)
occupation = np.zeros(size,dtype=complex)

#Define potentials

Vx = 0.5*(np.abs(x)-seperation/2.0)**2
Vk = (k**2)/2

#Define evolvers 
Uxh = np.exp(-1j*Vx*dt/2)
Ux = np.exp(-1j*Vx*dt)
Uf = np.exp(-1j*Vk*dt)

#define ground state of the H.O potentials
def Gauss(x): 
    return (np.pi**(-0.25))*np.exp((-x**2)/2.0) 

psi_l = Gauss(x+seperation/2)
psi_r = Gauss(x-seperation/2)

psi_0 = psi_l.astype(complex)

c = np.zeros(size,dtype=complex)
    
#evolution begins here

psi = psi_0
j = 0

#The operator we have to start off with

psi = psi*Uxh

psi_f = fft(psi)
psi_f = psi_f*Uf
psi = ifft(psi_f)

for i in range(nt):
        
    psi = psi*Ux
    
    psi_f = fft(psi)
    psi_f = psi_f*Uf
    psi = ifft(psi_f)
    
    if i % nmeasure == 0 :
        
        #Every time we measure , we have to finish with a Ux half time step
        
        psi_t = psi*Uxh
        
        c[j,0] =  np.abs(np.sum(np.conj(psi_t)*psi_l))*dx
        c[j,1] =  np.abs(np.sum(np.conj(psi_t)*psi_r))*dx
        
        j = j + 1
        
        print j
        
psi = psi*Uxh    
    
plt.plot(c*c)









