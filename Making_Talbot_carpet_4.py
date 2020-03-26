# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 21:49:13 2020

@author: Jay
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 18:56:49 2020

@author: Jay
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath

#this program is valid for a square grovved amplitude grating 


 
wavelength = 532*10**(-9) #wavelength of the laser
a = 20*10**(-6) #grating spacing
z_T = (a**2)/wavelength #Talbot length 
#h = wavelength #depth of phase modulation by the grating 
#rho = wavelength/(4*h*np.pi**2) #Bessel argument 
f = 0.1 

eta_limit = 4
zeta_limit = 8
eta_steps = 0.01
zeta_steps = 0.01
s = 0.25 #scaling factor

eta_grid = np.arange(0,eta_limit,eta_steps) #grid along direction of the grating
zeta_grid = np.arange(0,zeta_limit,zeta_steps) #grid along the propagation of the incident wave

number_of_orders = 10  #number of terms to add in 

size = (number_of_orders*2+1,np.size(eta_grid),np.size(zeta_grid)) #dimension of the wavefunction matrix
psi = np.zeros(size,dtype=complex)

for i in np.arange(-number_of_orders,number_of_orders+1,1): #storing values of wavefunction at each point for all orders
    if i == 0 :
        for j in range(np.size(eta_grid)):
            for k in range(np.size(zeta_grid)):
                psi[i,j,k] = f*complex(np.cos(np.pi*(2*i*eta_grid[j]-s*i*i*zeta_grid[k])),np.sin(np.pi*(2*i*eta_grid[j]-s*i*i*zeta_grid[k])))
        print i    
    else :     
        for j in range(np.size(eta_grid)):
            for k in range(np.size(zeta_grid)):
                psi[i,j,k] = (np.sin(np.pi*i*f)/(i*f))*complex(np.cos(np.pi*(2*i*eta_grid[j]-s*i*i*zeta_grid[k])),np.sin(np.pi*(2*i*eta_grid[j]-s*i*i*zeta_grid[k])))
        print i    

psi_summation = psi.sum(axis=0) #adding all orders to construct the final wavefunction

array_size = (np.size(eta_grid),np.size(zeta_grid)) #dimention of the intensity matrix

psi_absolute = np.zeros(array_size) #intensity matrix

for l in range(np.size(eta_grid)): #constructing intensity matrix
    for m in range(np.size(zeta_grid)):
        psi_absolute[l,m] = psi_summation[l,m]*np.conj(psi_summation[l,m])

fig, ax = plt.subplots() #the tick labelling is sensitive to the order of commands henceforth

plt.imshow(psi_absolute,cmap='jet') #create fig
plt.xticks(np.arange(0,zeta_limit/zeta_steps,400))  #plot initial ticks 
plt.yticks(np.arange(0,eta_limit/eta_steps,100))
labels_x = [item.get_text() for item in ax.get_xticklabels()] #prepare blanks for replacement ticks
labels_y = [item.get_text() for item in ax.get_yticklabels()]

for n in range(2): # fill up replacement tick arrays
    labels_x[n] = n 

for n in range(4):
    labels_y[n] = n
    
#labels[1] = 'Testing'
ax.set_xticklabels(labels_x) #set the replaced ticks
ax.set_yticklabels(labels_y)
#plt.savefig('talbot_carpet_amp_grating.svg')
plt.show()

