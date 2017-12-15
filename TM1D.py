# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 20:39:24 2017

@author: yhuang1
"""
#This program is used to calculate the statistical quantities in 1D random multilayer optic system based on transfer matrix
#The quantities includes average of logrithm of intensity <lnI(x)> and variance of logrithm of intensity var(ln(x))
from numpy import pi, conjugate, empty, linspace, concatenate, array, mean, log, var, floor, squeeze
from random import uniform
from cmath import sqrt,exp
from warnings import warn
from pylab import plot, xlabel, ylabel, title, xlim, ylim, show, gcf, legend
from matplotlib import rc
from time import clock
#transfer matrix calculated from wavenumber k1,k2 and layer thickness d1
#the 2-by-2 matrix is converted to 1-by-4 matrix
def TM(k1,k2,d1):
    phi=k1*d1
    phase=exp(-1j*phi)
    M=empty([1,4],complex)
    M[0,0]=(k2+k1)*phase
    M[0,1]=(k2-k1)*phase
    M[0,2]=conjugate(M[0,1])
    M[0,3]=conjugate(M[0,0])
    M=M*1/(2*k1)
    return M
# create the multilayer system 
def system_uniform_distribution(L_s,d_min,d_max,n1,n2): 
    # L_s is the length of the whole multilayer sample contained between the boundary, d_mean is the average value of layer thickness
    # sigma_d is the standard deviation of layer thickness, n_mean is the average refractive index of layers and sigma_n is the 
    # standard deviation of layer refractive index
    L=0
    d_system=[]
    n_system=[]
    #seed(1)
    i=0
    while(L<L_s):
        di=uniform(d_min,d_max)
        d_system.append(di)
        L+=di
        ni=n1+(n2-n1)*(i%2)
        n_system.append(ni)
        i+=1
    #cut the layer intersecting with the edge of end boundary
    d_system[-1]=di+L_s-L
    return d_system,n_system

def system_periodic(L_s,d1,d2,n1,n2):
    p=d1+d2
    if(p>L_s):
        warn('the size of the system is too small or the periodicity of the system is too large')
    L=0
    d_system=[]
    n_system=[]
    #seed(1)
    i=0
    while(L<L_s):
        ni=n1+(n2-n1)*(i%2)
        di=d1+(d2-d1)*(i%2)
        n_system.append(ni)
        d_system.append(di)
        L+=di
        i+=1
    #cut the layer intersecting with the edge of end boundary
    d_system[-1]=di+L_s-L
    return d_system,n_system
    
#create two boundaries of the sytem
def system_with_bdry(d_system,n_system,d_bdry1,d_bdry2,n_bdry1,n_bdry2):
    #d_system, n_system is the parameters of the original system and 
    #d_bdry1/2, n_bdry1/2 is the thickness and refractive index of each boundary layers
    d_system.insert(0,d_bdry1)
    d_system.append(d_bdry2)
    n_system.insert(0,n_bdry1)
    n_system.append(n_bdry2)
    return d_system,n_system

#find the serial number of layer where the each measuring point is
def layer_number_of_mp(x_interface,x_mp):
    delta=1E-12*x_interface[-1]
    if(max(x_mp)-max(x_interface)>delta or min(x_interface)-min(x_mp)>delta):
        warn('some measuring points are out of system')
    N_mp=len(x_mp)
    n_layer_x_mp=empty([N_mp,1],int)
    for n_mp in range(N_mp-1):
        n_interface=0
        while(x_interface[n_interface]-x_mp[n_mp]<=delta):
            n_interface+=1
        n_layer_x_mp[n_mp]=n_interface-1
    n_layer_x_mp[N_mp-1]=N_layer-1
    return n_layer_x_mp
#calculate out the wavevector at the propagating direction at a specific freqency
def freq2k(freq,n_system):
    k=[]
    a=22.86E-3 #width of waveguide
    vc=2.998E8 #velocity of light
    for i in range(len(n_system)):
        kx=pi*sqrt((2*n_system[i]*freq/vc)**2-(1/a)**2)
        k.append(kx)
    return array(k,complex)
## main program
t0=clock()
cv=2.998E8 #speed of light
L_s=2 #length of multilayer sample
n_s=50 #number of measuring points in the sample
d_bdry=[0.01,0.01]# thicknesses of two boundary layers
n_bdry=[1,1]# refractive indices of two boudary layers 
n_mp_b1=3#number of measuring points in left boundary
n_mp_b2=3#number of measuring points in right boudary
d_min=1E-5 #the minimum layer thickness 
d_max=0.06+1E-5 #the maximum layer thickness
d1,d2=0.0254,6.604E-3
n1=1    #refractive index of material1
n2=sqrt(6)    #refractive index of material2
# define measuring points
x_mp=concatenate((linspace(0,d_bdry[0],n_mp_b1),linspace(d_bdry[0]+L_s/(n_s+1),L_s+d_bdry[0]-L_s/(n_s+1),n_s),linspace(L_s+d_bdry[0],L_s+d_bdry[0]+d_bdry[1],n_mp_b2)),axis=0)
N_mp=len(x_mp)
N_config=int(1E4) #number of configurations
wl=1E-4 #wavelength
#freq=linspace(7,12,10)*1E9
freq=linspace(20,20.1,1)*1E9
N_freq=len(freq)
E=empty([N_config,N_freq,N_mp],complex) 
for n_config in range(N_config):
    #print(n_config)
    #create multilayer system with uniform distributon of layer thickness
    [d_system,n_system]=system_uniform_distribution(L_s,d_min,d_max,n1,n2)
    #[d_system,n_system]=system_periodic(L_s,d1,d2,n1,n2)
    [d_system,n_system]=system_with_bdry(d_system,n_system,d_bdry[0],d_bdry[1],n_bdry[0],n_bdry[1])
    
    #k=2*pi/wl*array(n_system)
    N_layer=len(d_system)
    #find the serial number of layer where the each measuring point is
    x_interface=[0]
    for i in range(N_layer):
        x_interface.append(x_interface[i]+d_system[i])
    x_interface=array(x_interface,float)
    n_layer_x_mp=layer_number_of_mp(x_interface,x_mp)
    
    for n_freq in range(N_freq):
        k=freq2k(freq[n_freq],n_system)
        #find the transfer matrix M[i,:] between every two adjacent layers i and i+1
        M=empty([N_layer-1,4],complex)
        for i in range(N_layer-1):
            M[i,:]=TM(k[i],k[i+1],d_system[i])
    
        #find the transfer matrix MN[i,:] between some layer i and the last layer N_layer-1
        MN=empty([N_layer,4],complex)
        MN[N_layer-1,:]=[1,0,0,1]
        for i in range(N_layer-2,-1,-1):
            MN[i,0]=M[i,0]*MN[i+1,0]+M[i,1]*MN[i+1,2]
            MN[i,1]=M[i,0]*MN[i+1,1]+M[i,1]*MN[i+1,3]
            MN[i,2]=M[i,2]*MN[i+1,0]+M[i,3]*MN[i+1,2]
            MN[i,3]=M[i,2]*MN[i+1,1]+M[i,3]*MN[i+1,3]
    
        #find the transfer matrix M0[i,:] between first layer 0 and some layer i
        M0=empty([N_layer,4],complex)    
        M0[0,:]=[1,0,0,1]
        for i in range(1,N_layer):
            M0[i,0]=M0[i-1,0]*M[i-1,0]+M0[i-1,1]*M[i-1,2]
            M0[i,1]=M0[i-1,0]*M[i-1,1]+M0[i-1,1]*M[i-1,3]
            M0[i,2]=M0[i-1,2]*M[i-1,0]+M0[i-1,3]*M[i-1,2]
            M0[i,3]=M0[i-1,2]*M[i-1,1]+M0[i-1,3]*M[i-1,3]
            #print(M0[i,:])
        

        #find the forward wave and backward wave for each layer 
        #based on formula (U[i],V[i])=(MN[i,0],MN[i,2])/(M0[i,0]*MN[i,0]+M0[i,1]*MN[i,2])
        U=empty([N_layer,1],complex)
        V=empty([N_layer,1],complex)
        for i in range(N_layer):
            cof=1/(M0[i,0]*MN[i,0]+M0[i,1]*MN[i,2])
            U[i]=cof*MN[i,0]
            V[i]=cof*MN[i,2]
    
        #find the total field at each measuring point
        #base on formula Ei(x)=U[i]*exp(1j*(x-xi))+V[i]*exp(-ij*(x-xi))
        for n_mp in range(N_mp):
            n_layer=n_layer_x_mp[n_mp]
            phase=exp(1j*k[n_layer]*(x_mp[n_mp]-x_interface[n_layer]))
            #print(k[n_layer],x_mp[n_mp],x_interface[n_layer],U[n_layer],V[n_layer],phase,U[n_layer]*phase+V[n_layer]*conjugate(phase))
            E[n_config,n_freq,n_mp]=U[n_layer]*phase+V[n_layer]*conjugate(phase)
    
t1=clock()
#the time used for caculation
print(t1-t0)

#from numpy import transpose;plot(freq,log(I[0,:,58]),freq,log(abs(transpose(E000))**2))

font = {'family' : 'Tahoma',
        'weight' : 'bold',
        'size'   : 25}
rc('font', **font)
#xlabel('freqency [Hz]')
#ylabel('ln(T)')
#title('transmission spectrum for periodic system')
#figure = gcf() # get current figure
#figure.set_size_inches(16, 12)
#show()
I=abs(E)**2
lnI=squeeze(log(I))
I=squeeze(I)
I_ave=mean(I,axis=0)
lnI_ave=mean(lnI,axis=0)
varlnI=var(lnI,axis=0)
plot(x_mp[n_mp_b1:-1-n_mp_b2],lnI_ave[n_mp_b1:-1-n_mp_b2],'.')
xlim((0,L_s))
ylim(ymax=0)
xlabel('x(measuring positons)')
ylabel('<lnI(x)>')
title('<lnI(x)> vs x')
figure = gcf() # get current figure
figure.set_size_inches(16, 12)
show()
plot(x_mp[n_mp_b1:-1-n_mp_b2],varlnI[n_mp_b1:-1-n_mp_b2],'.')
xlim((0,L_s))
xlabel('x(measuring positons)')
ylabel('var(lnI(x))')
title('var(lnI(x)) vs x')
figure = gcf() # get current figure
figure.set_size_inches(16, 12)
show()  
        
