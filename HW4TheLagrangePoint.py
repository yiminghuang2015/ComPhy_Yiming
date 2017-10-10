# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:22:56 2017

@author: yhuang1
"""

def f(r):
    # r is the distance between the Earth and satellite
    G=6.674E-11 #gravitational constant, unit:m^-3*Kg^-1*s^-2
    M=6.974E24 #the mass of the Earth, unit:Kg
    m=7.348E22 #the mass of the moon, unit:Kg
    R=3.844E8 #the distance between the Earth and the Moon, unit:m
    w=2.662E-6#the angular frequency of the Moon
    return G*(M/r**2-m/(R-r)**2)-w**2*r

r1=6.37E6
r2=3.84E8
d=1
if (f(r1)*f(r2)>0):
    print('can not find a solution between',r1,'and',r2)
    exit()
while (abs(r1-r2)>1E3):
    r=(r1+r2)/2
    f1=f(r1)
    f3=f(r)
    if(f1*f3>0):
        r1=r
    else:
        r2=r
print('The distance from the Earth to the Lagrange point is {:.3E}'.format(r),'m')
