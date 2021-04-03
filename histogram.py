import numpy as np
import h5py
from scipy.stats import skew, norm, kurtosis, lognorm
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick
import csv 
import scipy as sp
import math


######################## Input parameters #########################
files = 100
distributions = ['lognormal', 'uniform', 'normal']        #uniform, lognormal, normal



######################## Initialize processing files #########################
hist = []
fitos = []
par1 = []
par2 = []
min_v, max_v = -9, -5
x = np.logspace(min_v,max_v,1000)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
plt.xscale('log')

for counter,distribution in enumerate(distributions):
    
    hf = h5py.File(distribution+"/results_wu_"+distribution+"1.hdf5", 'r')
    Npoints = (np.array((hf.get('Coordinates')).get('Radius'))).size
    hf.close()
    Results = np.zeros(Npoints*files)
    ######################## Read files and join them in only one array #########################
    for j in range (0,files):
        filename = distribution+"/results_wu_"+distribution+str(j)+".hdf5"
        hf = h5py.File(filename, 'r')
        Results[j*Npoints:(j+1)*Npoints] = np.array((hf.get('Flow')).get('Results'))[:,5]
        hf.close()
        # print(j)
    mean = Results.mean()
    variance = np.var(Results)
    # par1.append(math.log((mean**2)/(math.sqrt(variance+mean**2))))
    # par2.append(math.log(variance/mean**2+1))
    fitos.append(sp.stats.lognorm.fit(Results,floc=0))
    print(Results.mean(), np.median(Results), np.var(Results))
    ######################## Histogram #########################
    plt.hist(Results, bins = x, label = distribution, density = True, histtype = 'step',linewidth = 3)
    # print(fitos[counter])
    print(fitos[counter][0],np.log(fitos[counter][2]))
pdf_fitted_automatic = sp.stats.lognorm.pdf(x, fitos[0][0], loc=fitos[0][1], scale=fitos[0][2])
plt.plot(x, pdf_fitted_automatic, '--', label = 'fit ln', color = 'blue')
pdf_fitted_automatic = sp.stats.lognorm.pdf(x, fitos[1][0], loc=fitos[1][1], scale=fitos[1][2])
plt.plot(x, pdf_fitted_automatic, '--', label = 'fit normal', color = 'orange')
pdf_fitted_automatic = sp.stats.lognorm.pdf(x, fitos[2][0], loc=fitos[2][1], scale=fitos[2][2])
plt.plot(x, pdf_fitted_automatic, '--', label = 'fit uniform', color = 'green')
# y = lognorm.pdf(x, fitos[0][0], fitos[0][1], fitos[0][2])
# plt.plot(x, y, label = 'hola', color = 'black',linewidth = 3)

sigma = fitos[0][0]
mu = math.log(fitos[0][2])
y = np.exp(-((np.log(x)-mu)**2/(2*sigma**2)))/(x*sigma*math.sqrt(2*math.pi))
plt.plot(x, y, label = 'chao', color = 'magenta',linewidth = 3)


# print(par1)
# print(par2)
plt.legend()
plt.show()


    


