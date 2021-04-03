from joblib import Parallel, delayed
from scipy import interpolate
import pandas as pd
import os, os.path
from wu import wu
import random
from variability_ranges import *
import h5py


######################## Input parameters #########################
files = 3
points = 100000
distribution = 'uniform'   #uniform, lognormal, normal
points_per_parameter = 1000

########### Be sure where to start #############
existing_files = 0
for file in os.listdir(distribution+'/.'):
    if file.startswith('results_wu'):
        existing_files += 1

########### Interpolate the state equation #############
properties = pd.read_csv('z,cg,n.csv')
p = properties['P'].values*1000000
t = properties['T'].values
f_z = interpolate.bisplrep(p, t, properties['Z'].values, s=2)
f_vis = interpolate.bisplrep(p, t, properties['Viscosity'].values/1000000, s=2)
f_cg = interpolate.interp1d(p, properties['cg'].values/1000000)

########### Create variability ranges #############
var_ranges = pd.read_csv('variability_ranges.csv').set_index('Variable')
var_lists = []
for item_par in var_ranges.index:
    var_lists.append(np.linspace(var_ranges.loc[item_par,'Min'], var_ranges.loc[item_par,'Max'], points_per_interval))

########### Create intervals and run #############
for j in range(existing_files,existing_files+files):
    if distribution == 'uniform':
        par_lists =  []
        for n_par, item_par in enumerate(var_ranges.index): par_lists.append(random.choices(var_lists[n_par], k =points))
    elif distribution == 'normal':
        pressure = np.random.normal(loc = 25250000, scale = 14289419, size = int(points*1.4))
        pressure = np.delete(pressure, np.where(pressure>50000000))
        pressure = np.delete(pressure, np.where(pressure<500000))

        radius = np.random.normal(loc = 5.1E-8, scale = 2.83E-8, size = int(points*1.4))
        radius = np.delete(radius, np.where(radius>1E-7))
        radius = np.delete(radius, np.where(radius<2E-9))

        porosity = np.random.normal(loc = 0.0525, scale = 0.0274, size = int(points*1.4))
        porosity = np.delete(porosity, np.where(porosity>0.1))
        porosity = np.delete(porosity, np.where(porosity<0.005))

        tortuosity = np.random.normal(loc = 4.3, scale = 0.866, size = int(points*1.4))
        tortuosity = np.delete(tortuosity, np.where(tortuosity>5.8))
        tortuosity = np.delete(tortuosity, np.where(tortuosity<2.8))

        op = np.random.normal(loc = 70500000, scale = 11258330, size = int(points*1.4))
        op = np.delete(op, np.where(op>90000000))
        op = np.delete(op, np.where(op<51000000))

        q = np.random.normal(loc = 0.035, scale = 0.012, size = int(points*1.4))
        q = np.delete(q, np.where(q>0.056))
        q = np.delete(q, np.where(q<0.014))

        t = np.random.normal(loc = 0.03, scale = 0.00577, size = int(points*1.4))
        t = np.delete(t, np.where(t>0.04))
        t = np.delete(t, np.where(t<0.02))

        k = np.random.normal(loc = 1.05, scale = 0.548, size = int(points*1.4))
        k = np.delete(k, np.where(k>2))
        k = np.delete(k, np.where(k<0.1))

        fd = np.random.normal(loc = 2.5, scale = 0.231, size = int(points*1.4))
        fd = np.delete(fd, np.where(fd>2.9))
        fd = np.delete(fd, np.where(fd<2.1))

        ih = np.random.normal(loc = 14000, scale = 1154.7, size = int(points*1.4))
        ih = np.delete(ih, np.where(ih>16000))
        ih = np.delete(ih, np.where(ih<12000))

        lp0 = np.random.normal(loc = 84500000, scale = 25114615, size = int(points*1.4))
        lp0 = np.delete(lp0, np.where(lp0>128000000))
        lp0 = np.delete(lp0, np.where(lp0<41000000))

        a0 = np.random.normal(loc = 1.19, scale = 0.1, size = int(points*1.4))
        a0 = np.delete(a0, np.where(a0>1.36))
        a0 = np.delete(a0, np.where(a0<1.02))

        a1 = np.random.normal(loc = 4, scale = 1.155, size = int(points*1.4))
        a1 = np.delete(a1, np.where(a1>6))
        a1 = np.delete(a1, np.where(a1<2))

        beta = np.random.normal(loc = 0.4, scale = 0.115, size = int(points*1.4))
        beta = np.delete(beta, np.where(beta>0.6))
        beta = np.delete(beta, np.where(beta<0.2))

        T = np.random.normal(loc = 405, scale = 39.26, size = int(points*1.4))
        T = np.delete(T, np.where(T>473))
        T = np.delete(T, np.where(T<337))

        # Create normal distributions of the parameters, these were created with mean and variance equal to those of the uniform distributions
        p_list = (pressure[:points]).tolist()           # Pore pressure                         Pa
        ro_list = (radius[:points]).tolist()            # Pore radius                           m
        por_list = (porosity[:points]).tolist()         # Porosity                              -
        tor_list = (tortuosity[:points]).tolist()       # Tortuosity                            -
        op_list = (op[:points]).tolist()                # Overburden pressure                   Pa
        q_list = (q[:points]).tolist()                  # q exponent porosity power law         -
        t_list = (t[:points]).tolist()                  # t exponent pore radius power law      -
        k_list = (k[:points]).tolist()                  # Blockage/Migration ratio              -
        fd_list = (fd[:points]).tolist()                # Fractal dimension of pore wall        -
        ih_list = (ih[:points]).tolist()                # Isosteric Heat                      J/mol
        lp0_list = (lp0[:points]).tolist()              # Langmuir pressure                     Pa
        a0_list = (a0[:points]).tolist()                # Alpha zero                            -
        a1_list = (a1[:points]).tolist()                # Alpha one                             -
        beta_list = (beta[:points]).tolist()            # Beta                                  -
        T_list = (T[:points]).tolist()                  # Temperature                           K
    elif distribution == 'lognormal':
        radius = np.random.lognormal(mean = -18.165, sigma = 1.01, size = int(points*1.4))
        radius = np.delete(radius, np.where(radius>1E-7))
        radius = np.delete(radius, np.where(radius<2E-9))
        ro_list = (radius[:points]).tolist()
        p_list = []
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
        beta_list = []
        T_list = []
        # Randomly sample the points
        for i in range(0,points):
            p_list.append(round(random.choice(p_range),1))
            por_list.append(random.choice(por_range))
            tor_list.append(random.choice(tor_range))
            op_list.append(random.choice(op_range))
            q_list.append(random.choice(q_range))
            t_list.append(random.choice(t_range))
            k_list.append(random.choice(k_range))
            fd_list.append(random.choice(fd_range))
            ih_list.append(random.choice(ih_range))
            lp0_list.append(random.choice(lp0_range))
            a0_list.append(random.choice(a0_range))
            a1_list.append(random.choice(a1_range))
            beta_list.append(random.choice(beta_range))
            T_list.append(random.choice(T_range))

    print('Im_in '+str(j))
    results = Parallel(n_jobs=-1)(delayed(wu)(p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,beta,T,f_z,f_vis,f_cg) for p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,beta,T in zip(par_lists[0],par_lists[1],par_lists[2],par_lists[3],par_lists[4],par_lists[5],par_lists[6],par_lists[7],par_lists[8],par_lists[9],par_lists[10],par_lists[11],par_lists[12],par_lists[13],par_lists[14]))
    print('I finished '+str(j))
    filename =distribution+"/results_wu_"+distribution+str(j)+".hdf5"
    hf = h5py.File(filename, "w")

    # Save each one of the files
    G0 = hf.create_group("Coordinates")
    for n_par, item_par in enumerate(var_ranges.index): G0.create_dataset(item_par, data=par_lists[n_par])
    G0.create_dataset("Z", data=np.array(results)[:,6])
    G0.create_dataset("cg", data=np.array(results)[:,7])
    G0.create_dataset("viscosity", data=np.array(results)[:,8])
    G1 = hf.create_group("Flow")
    G1.create_dataset("Results", data=results)
    hf.close()