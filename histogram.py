import numpy as np
import scipy as sp
import h5py
from scipy.stats import lognorm
from matplotlib import pyplot as plt 

######################## Input parameters #########################
files = 3
distributions = ['lognormal'] #, 'uniform', 'normal']        #uniform, lognormal, normal
min_v, max_v = -9, -5

######################## Initialize processing files #########################
fitos = []
x = np.logspace(min_v,max_v,1000)
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
plt.xscale('log')
for counter,distribution in enumerate(distributions):
    hf = h5py.File(distribution+"/results_wu_"+distribution+"1.hdf5", 'r')
    Npoints = (np.array((hf.get('Coordinates')).get('ro'))).size
    hf.close()
    Results = np.zeros(Npoints*files)
    ######################## Read files and join them in only one array #########################
    for j in range (0,files):
        filename = distribution+"/results_wu_"+distribution+str(j)+".hdf5"
        hf = h5py.File(filename, 'r')
        Results[j*Npoints:(j+1)*Npoints] = np.array((hf.get('Flow')).get('Results'))[:,5]
        hf.close()
    ######################## Fit a log-normal distribution, plot histogram and best fit #########################
    fitos.append(sp.stats.lognorm.fit(Results,floc=0))
    plt.hist(Results, bins = x, label = distribution, density = True, histtype = 'step',linewidth = 3)
    pdf_fitted_automatic = sp.stats.lognorm.pdf(x, fitos[0][0], loc=fitos[0][1], scale=fitos[0][2])
    plt.plot(x, pdf_fitted_automatic, '--', label = 'fit '+distribution)
plt.legend()
plt.show()


    


