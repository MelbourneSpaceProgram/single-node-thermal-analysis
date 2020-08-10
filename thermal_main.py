# -*- coding: utf-8 -*-
"""
******************Basic Thermal Analysis for a small cube-sat************************
----------------Melbourne Space Program, Melbourne, Australia-------------------------

Based on the works by: 
"Preliminary Thermal Analysis of Small Satellites" 
Casper Versteeg and David L. Cotten, Small Satellite Research Laboratory, University of Georgia, 

Author:
Raoul Mazumdar
PhD student in hypersonics propulsion,
Date: August 2020


Re     - Radius of Earth (m)
h      - Altitude (m)
Beta   - Orbital inclination angle (deg)
B_crit - Critical orbital inclination angle (deg)
a      - Albedo factor (-)
qir    - Heat flux due to IR
t      - Time(s)
sigma  - Stefan-Boltzmann constant
--------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams["font.family"] = "Times New Roman"

G = 6.67408e-11
Me = 5.972e24
sigma = 5.67e-8 # W.m^-2.K^-4
#
#
alpha = 0.96
epsilon = 0.90
#*************************** User input **********************************
h = 400*1000 # Alt
#Beta = 25    # Orbital inclination
Re = 6873*1000  # Radius of earth in meters
dt = 10   #  Simulation at dt time step (s)
#### Sat dimensions
Dia  = 21.11 / 100
Area = np.pi*(Dia**2)/4
Mc = 4
cp = 897 # Aluminium cp
#*************************************************************************
Beta_v  = [10,45,80]
colorkey =['k','b','g','r']
for vii in range(len(Beta_v)):
    Beta = Beta_v[vii]
    tp = np.sqrt(  (4*(np.pi**2)*((Re+h)**3)) / (G*Me)   )
    B_crit = np.arcsin(Re/(h+Re)) *(180/np.pi)
    
    if Beta < B_crit:
        nom = np.sqrt((h**2) + 2*Re*h)
        denom = (Re + h)*np.cos(Beta*(np.pi/180)) 
        fe = (1/180) * (np.arccos(nom/denom)*180/np.pi)
    else:
        fe = 0
        
        
    if Beta < 30:
        a = 0.14
        #
        qir = 228
        #
    elif Beta >= 30:
        a = 0.19
        #
        qir = 218
    
    #********* initialization model ************   
    def st_calc(x,tp,fe):
        t = x % tp
        if t < (tp/2)*(1-fe) or t > (tp/2)*(1+fe):
            st = 1
        else:
            st = 0
        return(st)   
        
    # Run the simulation for 1-year 
    
    tg_max  = 24*365*3600 # Seconds
    tg_max  = int(tp*10)
    
    qgen = 15.711 # Watts based on internal sat ssytems
    qsun = 1414 # W.m^-2 (Hot)
    #qsun = 1322 # w.m^-2 (cold)
    qsunv = [1414,1322]
    
    
    linekey = ['-',':']
    labelkey   = ['hot-solar value','cold-solar value']
    npts  =  int(tg_max / dt)
    
    
    for i in range(2):
        qsun = qsunv[i]
        time = []
        temp = []
        ti    = 0
        Tn = 293  #  Initial satellite temp
        for tx in range(npts):
            Qi  = qir*Area + (1+a)*qsun*st_calc( ti, tp, fe)*Area*alpha + qgen - (np.pi*(Dia**2)) *sigma*epsilon*(Tn**4) 
            Tn = Tn + Qi*dt/(cp*Mc) 
            #
            temp.append(Tn)
            time.append(ti)
            print(tx,'count',st_calc( ti, tp, fe))    
            ti = ti + dt
        y  = np.array(temp)-273.15
        x =  np.array(time)/tp
        plt.plot( x ,y, color=colorkey[vii], linestyle = linekey[i], label= '%s deg w/t %s'%(str(Beta),labelkey[i]) ) 

plt.xlabel('Time in orbit (Earth orbits)')
plt.ylabel('Temperature (Celcius)')
plt.legend(frameon=False,bbox_to_anchor=(1.03, 0.8))
plt.tight_layout()
plt.savefig('output_range.png', dpi=600)
plt.show()

    