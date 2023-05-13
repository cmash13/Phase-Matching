# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:16:16 2023

@author: Carter
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c



def extraordinary_index(wv,T):
#wvl and Temp can be either single valued or an array of equally spaced values (at least one must be a single value)
#Convention is for lambda to be in units of microns in the ne equation!!!!!!!!!
#Based off this paper DOI: 10.1007/s00340-008-2998-2
    T0 = 24.5
  
    a1 = 5.756
    b1 = 2.860*10**(-6)

    a2 = 0.0983
    b2 = 4.700*10**(-8)

    a3 = 0.2020
    b3 = 6.113*10**(-8)

    a4 = 189.32
    b4 = 1.516*10**(-4)

    a5 = 12.52
    b5 = 9.62119*10**(-10)

    a6 =  1.32*10**(-2)
    f = (T - 24.5)*(T + 570.82)
    
    ne_t1 = a1
    ne_t2 = b1 * f
    ne_t3 = (a2 + (b2*f))/(wv**2 - (a3 + (b3*f))**2)
    ne_t4 = (a4 + (b4*f))/(wv**2 - a5**2)
    ne_t5 = a6*(wv**2)
    #ne_t5 = a6*wv**2
    ne = np.sqrt((ne_t1) + (ne_t2) + (ne_t3) + (ne_t4) - (ne_t5))
    
    
    return ne, wv, T

def extraordinary_index_omega(omega,T):
    
    T0 = 24.5
  
    a1 = 5.756
    b1 = 2.860*10**(-6)

    a2 = 0.0983
    b2 = 4.700*10**(-8)

    a3 = 0.2020
    b3 = 6.113*10**(-8)

    a4 = 189.32
    b4 = 1.516*10**(-4)

    a5 = 12.52
    b5 = 9.62119*10**(-10)

    a6 =  1.32*10**(-2)
    f = (T - 24.5)*(T + 570.82)
    
    wv = (2*np.pi*c)/(omega)
    
    ne_t1 = a1
    ne_t2 = b1 * f
    ne_t3 = (a2 + (b2*f))/(wv**2 - (a3 + (b3*f))**2)
    ne_t4 = (a4 + (b4*f))/(wv**2 - a5**2)
    ne_t5 = a6*(wv**2)
    #ne_t5 = a6*wv**2
    ne = np.sqrt((ne_t1) + (ne_t2) + (ne_t3) + (ne_t4) - (ne_t5))
    
    return ne,wv
 
    
def dndw(wv,T):
    #math function for helping to define group velocity.... this is taking care of the derivative of ne WRT omega
    omega = (2*np.pi*c)/(wv)
    
    T0 = 24.5

    a2 = 0.0983
    b2 = 4.700*10**(-8)

    a3 = 0.2020
    b3 = 6.113*10**(-8)

    a4 = 189.32
    b4 = 1.516*10**(-4)

    a5 = 12.52
    b5 = 9.62119*10**(-10)

    a6 =  1.32*10**(-2)
    f = (T - 24.5)*(T + 570.82)
    
    dndw_t_all = (8*(np.pi**2)*(c**2))/(omega**3)
    dndw_t_common = 4*(c**2)*(np.pi**2)/(omega**2)
    
    dndw_t1 = a6
    dndw_t2 = (a4 + b4*f)/((dndw_t_common - a5**2)**2)
    dndw_t3 = (a2 + b2*f)/((dndw_t_common - (a3 + b3*f)**2)**2)
    
    dndw = dndw_t_all*(dndw_t1 + dndw_t2 + dndw_t3)
    
    return dndw.astype(float)
 
def Vg_extraordinary(wv,T):
    #calculating group velcoity associated with extraordinary index of refraction
    
    omega = (2*np.pi*c)/(wv)
    
    dkdw_t1,junk = extraordinary_index_omega(omega, T)
    dkdw_t2 = (dndw(wv,T)*omega)*(1/c)
    
    dkdw = (dkdw_t1/c) + dkdw_t2
    
    vg = 1/dkdw
    
    return vg
   
    

def Temp_Poling(polingPeriod,T):
    
        #temperature should be in Kelvin!!!
    
        alpha = 1.54E-5
        beta = 5.3E-9
        PolingOfT = polingPeriod*(1 + alpha*(T-298) + beta*((T-298)**2))
        
        return PolingOfT,T


def deltaK_TypeZeroSFG(n1,wv1,n2,wv2,n3,wv3,poling):
    #poling is defind in microns already
    #does not include temperature dependece of poling period
    #n1,wv1 is SFG field that is desired
    
    
    
    dK_t1 = n1/wv1
    dK_t2 = n2/wv2
    dK_t3 = n3/wv3
    dK_t4 = 1/(poling)
    
    dK = 2*np.pi*(dK_t1 - dK_t2 - dK_t3 - dK_t4)
    return dK,poling
    
#---------------------------------------------------------------------    
"plotting wavevector mismatch deltaK as a function of poling period for frequency doubling via SFG from 1-2 micron"   

def deltaKvsPolingSHG(wavelength_in,PP,T):
   

    for i in wavelength_in:
        ne_in,wv_in,T = extraordinary_index(i,T)
        ne_out,wv_out,T = extraordinary_index(i/2,T)
        dK,poling = deltaK_TypeZeroSFG(ne_out,wv_out,ne_in,wv_in,ne_in,wv_in,PP)
        plt.plot(poling,dK,label= '%.2f \u03BCm'% wv_in)
        plt.axhline(y=0,linestyle='dashed')
        plt.legend()
        plt.xlabel("Poling Period (\u03BCm)")
        plt.ylabel("\u0394K")


#-----------------------------------------------------------------------------
def GroupVelocityMismatch(wv_in,T):
    
    vg_in = Vg_extraordinary(wv_in, T)
    vg_SHG = Vg_extraordinary(wv_in/2, T)
    
    gvmm = (1/vg_SHG) - (1/vg_in)

    return gvmm



#------------------------------------------------------------------------
def Intensity_SHG(gvmm,wv,dk,L,I_in):
    prefactor = 1
    I_SHG_t1 = (((2*np.pi*c)/wv)*gvmm)*(L/2)
    I_SHG_t2 = (dk*L)/2
    
    I_SHG = prefactor*((np.sinc(I_SHG_t1 - I_SHG_t2))**2)*(I_in**2)
    
    return I_SHG,wv,dk,L,I_in


#--------------------------------------------------------------------------------
wv_in = np.linspace(0.5,5,num=1000000)
T = 25
polingPeriod = 18.5
L_crystal = 1E-3
I_in = 1
ne_in,junk,junkT = extraordinary_index(wv_in, T) 
ne_SHG,junk,junk = extraordinary_index(wv_in/2, T)


gvmm = GroupVelocityMismatch(wv_in, T)
dK,junk = deltaK_TypeZeroSFG(ne_SHG,wv_in/2,ne_in,wv_in,ne_in,wv_in,polingPeriod)
I_SHG,wv,dk,L,I_in = Intensity_SHG(gvmm,(wv_in)*10**(-6),dK,L_crystal,I_in)

plt.plot(wv_in,I_SHG)
#deltaKvsPolingSHG(np.linspace(1,2,num=10),np.linspace(5,30,num=1000),200)

#ne,wv = extraordinary_index_omega(1, 25)
"""   
ne1,wv1,T = extraordinary_index(np.linspace(0.5,4), 20)

plt.plot(wv1,ne1)
plt.ylabel("Refractive Index")
plt.xlabel("Wavelength (\u03BCm)")
ax = plt.gca()
ax.get_yaxis().get_major_formatter().set_useOffset(False)
"""