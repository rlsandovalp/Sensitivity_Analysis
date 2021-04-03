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
distribution = 'normal'   #uniform, lognormal, normal
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
        par_lists = []
        for n_par, item_par in enumerate(var_ranges.index): 
            mean = (var_ranges.loc[item_par,'Min'] + var_ranges.loc[item_par,'Max'])/2
            std = (var_ranges.loc[item_par,'Max'] - var_ranges.loc[item_par,'Min'])/(12**0.5)
            parametro_gen = np.random.normal(loc = mean, scale = std, size = int(points*1.4))
            parametro_gen = np.delete(parametro_gen, np.where(parametro_gen>var_ranges.loc[item_par,'Max']))
            parametro_gen = np.delete(parametro_gen, np.where(parametro_gen<var_ranges.loc[item_par,'Min']))
            par_lists.append(parametro_gen[:points].tolist())

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