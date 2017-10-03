# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 17:21:29 2017

@author: yhuang1
"""
#5.11 Plane Wave
#import sys 
#sys.path.append(r'C:\Users\yhuang1\Dropbox\Courses\Computational Physics')
#import gaussxw

from numpy import ones,copy,cos,tan,pi,linspace

def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

import numpy as np
import pylab
wl=1  #wavelength(m)
z=3
x=np.linspace(-5,5,1000)
Ir=np.zeros(len(x))
t0,w0=gaussxw(50)
a=0
for i in range(len(x)):
    u=x[i]*np.sqrt(2/wl/z)
    b=u
    c,s=0,0
    for t,w in zip(t0,w0):
        tt=1/2*(b-a)*t+1/2*(b+a)
        ww=1/2*(b-a)*w
        c+=np.cos(1/2*np.pi*tt**2)*ww
        s+=np.sin(1/2*np.pi*tt**2)*ww
    Ir[i]=1/8*((2*c+1)**2+(2*s+1)**2)
pylab.xlabel('x(m)')
pylab.ylabel('I/I0')
pylab.plot(x,Ir)
pylab.show()