from joblib import Parallel, delayed
from scipy import interpolate
import numpy as np
import pandas as pd
import os, os.path
from wu import wu
import random
import h5py


######################## Input parameters #########################
files = 3
points = 100000
distribution = 'lognormal'   #uniform, lognormal, normal
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
for item_par in var_ranges.index: var_lists.append(np.linspace(var_ranges.loc[item_par,'Min'], var_ranges.loc[item_par,'Max'], points_per_parameter))

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
        par_lists = []
        mean = (var_ranges.loc['p','Min'] + var_ranges.loc['p','Max'])/2
        std = (var_ranges.loc['p','Max'] - var_ranges.loc['p','Min'])/(12**0.5)
        pressure = np.random.normal(loc = mean, scale = std, size = int(points*1.4))
        pressure = np.delete(pressure, np.where(pressure>var_ranges.loc['p','Max']))
        pressure = np.delete(pressure, np.where(pressure<var_ranges.loc['p','Min']))
        par_lists.append(pressure[:points].tolist())
        radius = np.random.lognormal(mean = -18.165, sigma = 1.01, size = int(points*1.4))   # This parameters correspond to a given mean specified in the article
        radius = np.delete(radius, np.where(radius>1E-7))
        radius = np.delete(radius, np.where(radius<2E-9))
        par_lists.append(radius[:points].tolist())
        a = var_ranges.index.tolist()
        a.remove('p')
        a.remove('ro')
        for n_par, item_par in enumerate(a): 
            mean = (var_ranges.loc[item_par,'Min'] + var_ranges.loc[item_par,'Max'])/2
            std = (var_ranges.loc[item_par,'Max'] - var_ranges.loc[item_par,'Min'])/(12**0.5)
            parametro_gen = np.random.normal(loc = mean, scale = std, size = int(points*1.4))
            parametro_gen = np.delete(parametro_gen, np.where(parametro_gen>var_ranges.loc[item_par,'Max']))
            parametro_gen = np.delete(parametro_gen, np.where(parametro_gen<var_ranges.loc[item_par,'Min']))
            par_lists.append(parametro_gen[:points].tolist())

    print('Im_in '+str(j))
    results = Parallel(n_jobs=-1)(delayed(wu)(p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,beta,T,f_z,f_vis,f_cg) for p,ro,por,tor,op,q,t,k,fd,ih,lp0,a0,a1,beta,T in zip(par_lists[0],par_lists[1],par_lists[2],par_lists[3],par_lists[4],par_lists[5],par_lists[6],par_lists[7],par_lists[8],par_lists[9],par_lists[10],par_lists[11],par_lists[12],par_lists[13],par_lists[14]))

    # Save each one of the files
    filename =distribution+"/results_wu_"+distribution+str(j)+".hdf5"
    hf = h5py.File(filename, "w")
    G0 = hf.create_group("Coordinates")
    for n_par, item_par in enumerate(var_ranges.index): G0.create_dataset(item_par, data=par_lists[n_par])
    G0.create_dataset("Z", data=np.array(results)[:,6])
    G0.create_dataset("cg", data=np.array(results)[:,7])
    G0.create_dataset("viscosity", data=np.array(results)[:,8])
    G1 = hf.create_group("Flow")
    G1.create_dataset("Results", data=results)
    hf.close()