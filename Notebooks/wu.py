import math
import pandas as pd
from scipy import interpolate

athm_pr = 100000                           # [Pa]
dm = 3.8E-10                               # [-]
R = 8.314462                               # [J/(mol - K)]
M = 1.6E-2                                 # [kg/mol]
Na = 6.0221415E23                          # [1/mol]
b = -1
dpdl = 0.1                                 # [MPa/m]

def wu(p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,beta,T,f_z,f_vis,f_cg,f_rho):
    Z = interpolate.bisplev(p,T,f_z)
    vis = interpolate.bisplev(p,T,f_vis)
    rho = interpolate.bisplev(p,T,f_rho)
    cg = f_cg(p)
    dzdp = (interpolate.bisplev(p+10000,T,f_z)-interpolate.bisplev(p-10000,T,f_z))/20000
    
    porosity = por*((op-p)/athm_pr)**-q  
    radius = ro*((op-p)/athm_pr)**-t 
    kn = (vis/(p*2*radius))*(math.pi*Z*R*T/(2*M))**0.5
    wv = 1/(1+kn)
    wk = 1-wv
    lp = lp0 * math.exp(-ih/(R*T))
    theta = (p/Z)/(lp+p/Z)
    rad = theta*dm    ####AAAAAAAAAAAAAAAAAAAAAAAAA
    fc1 = (porosity/tor)*(1-rad/radius)**2
    fc2 = fc1 * (((1-rad/radius)**-2)-1)
    alpha = a0*2*math.atan(a1*kn**beta)/math.pi 
    B1 = wv*fc1*radius**2*p*M*(1+alpha*kn)*(1+4*kn/(1-b*kn))/(8*vis*Z*R*T)
    B2 = wk*fc1*2*radius*(dm/radius)**(fd-2)*((8*Z*M/(math.pi*R*T))**0.5)*p*cg/(3*Z)
    if k>=1:
        gam = 0
    else:
        gam = 1 
    ds0 = 8.29E-7 * T**0.5 * math.exp(-ih**0.8/(R*T)) # m2/s
    ds = ds0 * ((1-theta) + 0.5*k*theta*(2-theta) + gam*0.5*k*(1-k)*(theta**2))/((1-theta+0.5*k*theta)**2) # m2/s
    Csc = 4*theta*M/(math.pi*dm**3*Na) # kg*mol/mol*m^3    kg/m^3
    B3 = fc2*ds*Csc/(p) # m^2*kg*s^2*m^2/s*m^3*kg*m = s
    D1 = (B1*R/M)/(1/(T*Z)+p*dzdp/(T*Z**2))
    D2 = (B2*R/M)/(1/(T*Z)+p*dzdp/(T*Z**2))
    D3 = (B3*R/M)/(1/(T*Z)+p*dzdp/(T*Z**2))
    DT = D1+D2+D3
    J = (B1+B2+B3)*dpdl*86400*365/1000
    return (B1, B2, B3, J, D1, D2, D3, DT, kn, dzdp, Z, cg, vis, rho)