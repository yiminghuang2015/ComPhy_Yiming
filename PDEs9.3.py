# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 21:45:18 2017

@author: Yiming
"""

import numpy as np
import pylab as pl

#constants
M=100
v=1
target=1E-6
#initial value of potential and boundary conditions
phi=np.zeros([M+1,M+1],float)
phi[20:81,20]=v
phi[20:81,80]=-v
phiprime=np.empty([M+1,M+1],float)
#main loop
delta=1
while delta>target:
    for i in range(M+1):
        for j in range(M+1):
            if i==0 or i==M or j==0 or j==M or (j==20 and i>=20 and i<=80) or (j==80 and i>=20 and i<=80):
                phiprime[i,j]=phi[i,j]
            else:
                phiprime[i,j]=(phi[i+1,j]+phi[i-1,j]+phi[i,j+1]+phi[i,j-1])/4
                #find the maximum difference from old values
    delta=np.amax(abs(phi-phiprime))
    phi,phiprime=phiprime,phi
#make a plot
pl.imshow(phi)
pl.gray()
pl.show()

        