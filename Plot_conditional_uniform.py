import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


distribution = 'lognormal'
sm = 'Variance'

if sm == 'Mean':
    file_sm = 'CM'
    sm_pos = 0
    range_plot = (1E-9,1E-6)
elif sm == 'Variance':
    file_sm = 'CV'
    sm_pos = 1
    range_plot = (1E-18,1E-13)
elif sm == 'Skewness':
    file_sm = 'CS'
    sm_pos = 2
    range_plot = (1E-9,1E-6)
elif sm == 'Kurtosis':
    file_sm = 'CK'
    sm_pos = 3
    range_plot = (1E-9,1E-6)

unc = np.loadtxt('Results/'+distribution+'_SM.txt')[sm_pos]*np.ones(100)
fig = plt.figure(figsize=(10,6))
axes = fig.add_subplot(111)
size_label = 18

CMoment = np.loadtxt('Results/'+distribution+'_'+file_sm+'.txt')
axes.plot(unc, linewidth = 3, c = 'gray', label = sm)
axes.plot(CMoment[0,:], '.', c = 'black', label = r'$r_o$')
axes.plot(CMoment[1,:], '.', c = 'purple', label = r'$\phi_o$')
axes.plot(CMoment[2,:], '.', c = 'red', label = r'$p$')
axes.plot(CMoment[3,:], '.', c = 'blue', label = r'$\tau$')
axes.plot(CMoment[4,:], '.', c = 'green', label = r'$T$')
axes.set_ylim(range_plot)
axes.set_yscale('log')
axes.set_xlabel('Normalized parameter', fontsize = size_label)
axes.set_ylabel('Conditional '+sm, fontsize = size_label)
axes.tick_params(axis = "y", which = 'both', right = True, direction = "in", labelsize = size_label)
axes.xaxis.set_minor_locator(AutoMinorLocator())
axes.tick_params(axis = "x", which = 'both', top = True, direction = "in", labelsize = size_label)
axes.legend(prop={'size': 15})

plt.show()