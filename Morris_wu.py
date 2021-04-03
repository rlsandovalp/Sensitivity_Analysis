import numpy as np
import pandas as pd
import random
import math
import h5py

# Select the number of random points (r)
points = 1000

# Wu model
def wu(p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,betha):
    Z = (properties.loc[properties['T']==T][properties['P']==round(p/1000000,1)]).iloc[0]['Z']                          # -
    cg = (properties.loc[properties['T']==T][properties['P']==round(p/1000000,1)]).iloc[0]['cg']/1000000                # 1/Pa
    vis = (properties.loc[properties['T']==T][properties['P']==round(p/1000000,1)]).iloc[0]['Viscosity']/1000000        # Pa-s
    porosity = por*((op-p)/athm_pr)**-q
    radius = ro*((op-p)/athm_pr)**-t
    kn = (vis/(p*2*radius))*(math.pi*Z*R*T/(2*M))**0.5
    wv = 1/(1+kn)
    wk = 1/(1+1/kn)
    lp = lp0 * math.exp(-ih/(R*T))
    theta = (p/Z)/(lp+p*Z)
    rad = theta*dm
    fc1 = (porosity/tor)*(1-rad/radius)**2
    fc2 = fc1 * (((1-rad/radius)**-2)-1)
    alpha = a0*2*math.atan(a1*kn**betha)/math.pi 
    JJv = fc1*radius**2*p*(1+alpha*kn)*(1+4*kn/(1-b*kn))/(8*vis)
    JJk = fc1*2*radius*(dm/radius)**(fd-2)*((8*Z*R*T/(math.pi*M))**0.5)*p*cg/3
    if k>=1:
        gam = 0
    else:
        gam = 1
    ds0 = 8.29E-7 * T**0.5 * math.exp(-ih**0.8/(R*T)) 
    ds = ds0 * ((1-theta) + 0.5*k*theta*(2-theta) + gam*0.5*k*(1-k)*(theta**2))/((1-theta+0.5*k*theta)**2)
    Csc = 4*theta*M/(math.pi*dm**3*Na)
    JJs = fc2*ds*Csc*R*T*Z/(p*M)                 
    dj1 = (wv*JJv + wk*JJk + JJs)*M*0.1/(R*T*Z)*86400*365/1000
    return dj1

# Define constants
global athm_pr,dm,R,M,Na,b,a,T,lp0,a0,a1,betha,properties
athm_pr = 100000                        # MPa
dm = 3.8E-10                            # meters
R = 8.314462                            # J/(mol - K)
M = 1.6E-2                              # kg/mol
Na = 6.0221415E23                       # 1/mol
b = -1
a = 0
T = 435                                 # K
properties = pd.read_csv('z,cg,n.csv')

# Allocate lists
global p_list,ro_list,por_list,tor_list,op_list,q_list,t_list,k_list,fd_list,ih_list
p_list = []
ro_list = []
por_list = []
tor_list = []
op_list = []
q_list = []
t_list = []
k_list = []
fd_list = []
ih_list = []
lp0_list = []
a0_list = []
a1_list = []
betha_list = []

# Define ranges of variability of the parameters
p_range, dp = np.linspace (500000, 50000000, 100, retstep=True)           # Pore pressure                         Pa
ro_range, dro = np.linspace(2E-9, 1E-8, 100, retstep=True)                # Pore radius                           m
por_range, dpor = np.linspace(0.005, 0.1, 100, retstep=True)              # Porosity                              -
tor_range, dtor = np.linspace(2.8, 5.8, 100, retstep=True)                # Tortuosity                            -
op_range, dop = np.linspace(51000000, 90000000, 100, retstep=True)        # Overburden pressure                   Pa
q_range, dq = np.linspace(0.014, 0.056, 100, retstep=True)                # q exponent porosity power law         -
t_range, dt = np.linspace(0.02, 0.04, 100, retstep=True)                  # t exponent pore radius power law      -
k_range, dk = np.linspace(0.1, 2, 100, retstep=True)                      # Blockage/Migration ratio              -
fd_range, dfd = np.linspace(2.1, 2.9, 100, retstep=True)                  # Fractal dimension of pore wall        -
ih_range, dih = np.linspace(12000, 16000, 100, retstep=True)              # Isosteric heat of the geomaterial    J/mol-K
lp0_range = np.linspace(41000000, 128000000, 1000)                                      # Langmuir pressure                     Pa
a0_range = np.linspace(1.02, 1.36, 1000)                                    # Alpha zero                            -
a1_range = np.linspace(2, 6, 1000)                                          # Alpha one                             -
betha_range = np.linspace(0.2, 0.6, 1000)                                    # Beta                                  -


# Randomly sample the points
for i in range(0,points):
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    p_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    ro_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    por_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    tor_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    op_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    q_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    t_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    k_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    fd_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    ih_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    lp0_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    a0_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    a1_list.append(algo)
    algo = []
    algo.append(random.randrange(98)+1)
    algo.append(algo[0]+1)
    algo.append(algo[0]-1)
    betha_list.append(algo)


e_effects = np.zeros((points,14,2))
basis = np.zeros(points)

for i in range (0,points):
    basis[i] = wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])
    # Pore pressure
    e_effects[i,0,0] = (wu(p_range[p_list[i][1]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,0,1] = (basis[i]-wu(p_range[p_list[i][2]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # Pore radius
    e_effects[i,1,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][1]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,1,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][2]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # Porosity
    e_effects[i,2,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][1]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,2,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][2]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # Tortuosity
    e_effects[i,3,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][1]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,3,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][2]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # OP
    e_effects[i,4,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][1]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,4,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][2]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # q
    e_effects[i,5,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][1]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,5,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][2]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # t
    e_effects[i,6,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][1]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,6,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][2]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # k
    e_effects[i,7,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][1]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,7,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][2]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # fd
    e_effects[i,8,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][1]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,8,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][2]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # ih
    e_effects[i,9,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][1]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,9,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][2]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # lp0
    e_effects[i,10,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][1]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,10,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][2]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # a0
    e_effects[i,11,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][1]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,11,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][2]], a1_range[a1_list[i][0]], betha_range[betha_list[i][0]]))/0.01
    # a1
    e_effects[i,12,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][1]], betha_range[betha_list[i][0]])-basis[i])/0.01
    e_effects[i,12,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][2]], betha_range[betha_list[i][0]]))/0.01
    # betha
    e_effects[i,13,0] = (wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][1]])-basis[i])/0.01
    e_effects[i,13,1] = (basis[i]-wu(p_range[p_list[i][0]], ro_range[ro_list[i][0]], por_range[por_list[i][0]], tor_range[tor_list[i][0]], op_range[op_list[i][0]], q_range[q_list[i][0]], t_range[t_list[i][0]], k_range[k_list[i][0]], fd_range[fd_list[i][0]], ih_range[ih_list[i][0]], lp0_range[lp0_list[i][0]], a0_range[a0_list[i][0]], a1_range[a1_list[i][0]], betha_range[betha_list[i][2]]))/0.01


morr_ind_avr = np.zeros(14)
morr_ind_var = np.zeros(14)

for i in range(0,14):
    morr_ind_avr[i] = np.mean(np.abs(e_effects[:,i,:]))
    morr_ind_var[i] = np.var(np.abs(e_effects[:,i,:]))

print(morr_ind_avr)
print(morr_ind_var)