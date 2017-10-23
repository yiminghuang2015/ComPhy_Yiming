# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 16:39:52 2017

@author: Administrator
"""

import numpy as np
import pylab as pl
#global variables for t, x, y
globalt=[0]
globalx=[0]
globaly=[0]
def f(r,t):
    a,b=1,3
    x=r[0]
    y=r[1]
    fx=1-(b+1)*x+a*x**2*y
    fy=b*x-a*x**2*y
    return np.array([fx,fy],float)
def step(r,t,H):
    delta=1E-10 #accuracy per unit time for both x and y
    n=1
    r1=r+0.5*H*f(r,t)
    r2=r+H*f(r1,t+0.5*H)
    R1=np.empty([1,2],float)
    R1[0]=0.5*(r1+r2+0.5*H*f(r2,t+H))
    error=2*H*delta
    while (error>H*delta and n<=8):
        n+=1
        h=H/n
        #Modified midpoint method
        r1=r+0.5*h*f(r,t)
        r2=r+h*f(r1,t+0.5*H)
        
        for i in range(n-1):
            r1+=h*f(r2,t+(i+1)*h)
            r2+=h*f(r1,t+(i+1.5)*h)
            
        #R2 is the perivous row of extrapolation estimates
        # R1 is the latest row of extrapolation estimates
        R2=R1
        
    
        R1=np.empty([n,2],float)
        R1[0]=0.5*(r1+r2+0.5*h*f(r2,t+H))
        for m in range(1,n):
            epsilon=(R1[m-1]-R2[m-1])/((n/(n-1))**(2*m)-1)
            R1[m]=R1[m-1]+epsilon
        error=max(abs(epsilon))
        # in case that the error or fx or fy is too large to overflow
        if (error>1E50):
            break
    rr=R1[n-1]
    
    # recursive method 
    if (error>H*delta):
        r11=step(r,t,H/2)
        r22=step(r11,t+H/2,H/2)
        return np.array([r22[0],r22[1]],float)
    else:
        globalt.append(t+H)
        globalx.append(rr[0])
        globaly.append(rr[1])
        #pl.plot([t,t],[rr[0],rr[1]],'o')
        return np.array([rr[0],rr[1]],float)

step([0,0],0,20)
#make a plot
plot1=pl.plot(globalt,globalx,'o--',label='$x$')
plot2=pl.plot(globalt,globaly,'*--',label='y')
pl.xlabel("time")
pl.ylabel('x(t) and y(t)')
pl.legend()
#pl.savefig('C:\\Users\\yhuang1\\Documents\\Spyder\\foo.jpg', bbox_inches='tight')
pl.show()
     
     
     
     
     
     


