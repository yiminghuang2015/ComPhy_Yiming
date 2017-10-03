# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#5.9 Heat capacity of a solid
#import sys 
#sys.path.append(r'C:\Users\yhuang1\Dropbox\Courses\Computational Physics')
#import gaussxw
import numpy as np
import pylab
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

def cv(T):
    if T<=0:
        print('Error! The temperature should be higher than 0')
        return None
    N=50 # number of sample points used in Gaussian quadrature
    n_d=6.022E28 # number density of solid aluminum in unit m^-3
    t_D=428 # Debye temperature in unit K
    V=1000E-6 #solid volume in unit m^-3
    kb=1.3806E-23 # Boltzmannn constant in unit m^2*Kg*s^-2*K^-1
    c_t=9*V*n_d*kb/(t_D**3)*np.exp(4) # total constant for Cv
    t,w=gaussxwab(N,0,t_D/T)
    I=np.empty(1,float)
    for ti,wi in zip(t,w):
        I+=wi*ti**4/(np.exp(ti)-1)**2
    Cv=I*c_t
    return Cv


T1=5  #minimum temperature
T2=500 #maximum temperature
n=1000 #number of points in range [T1,T2]
T=np.linspace(T1,T2,n)
Cv=np.empty(n,float)
for i in range(len(T)):
    Cv[i]=cv(T[i])
pylab.xlabel('Temperature T (K)')
pylab.ylabel('Heat capacity Cv for 1000 cm^3 aluminum (J*K^-1)')
pylab.plot(T,Cv)
pylab.show()
    
