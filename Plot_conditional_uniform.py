import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


distribution = 'Uniform'

sm_to_plot = ['mean', 'variance', 'skew', 'kurt']
variables_to_plot = ['Unconditioned', 'Radius', 'Porosity', 'Pressure', 'Tortuosity', 'T']

CM = np.loadtxt('Results/'+distribution+'_CM.txt')
CV = np.loadtxt('Results/'+distribution+'_CV.txt')



fig, axes = plt.subplots(nrows = 2, ncols = 2, dpi = 250)
size_points = 1
size_label = 7
constant = 3


### MEAN
axes[0,0].text(0.9, 0.85, '(a)', transform=axes[0,0].transAxes, fontsize = size_label+constant, weight='bold')
axes[0,0].plot(CM[0,:], c = 'black', label = 'Radius')
axes[0,0].plot(CM[1,:], c = 'purple', label = 'Porosity')
axes[0,0].plot(CM[2,:], c = 'red', label = 'Pressure')
axes[0,0].plot(CM[3,:], c = 'blue', label = 'Tortuosity')
axes[0,0].plot(CM[4,:], c = 'green', label = 'Temperature')
# axes[0,0].plot(points_x,mean.loc['Unconditioned'].values, linewidth = size_points*2, c = 'gray')
axes[0,0].set_ylim(1E-9,1E-6)
axes[0,0].set_yscale('log')
axes[0,0].set_ylabel('Conditional Mean', fontsize = size_label)
axes[0,0].tick_params(axis = "y", which = 'both', right = True, direction = "in", labelsize = size_label)
axes[0,0].xaxis.set_minor_locator(AutoMinorLocator())
axes[0,0].tick_params(axis = "x", which = 'both', top = True, direction = "in", labelsize = size_label)
axes[0,0].legend()

### VARIANCE
axes[0,1].text(0.9, 0.85, '(b)', transform=axes[0,1].transAxes, fontsize = size_label+constant, weight='bold')
axes[0,1].plot(CV[0,:], c = 'black')
axes[0,1].plot(CV[1,:], c = 'purple')
axes[0,1].plot(CV[2,:], c = 'red')
axes[0,1].plot(CV[3,:], c = 'blue')
axes[0,1].plot(CV[4,:], c = 'green')
# axes[0,1].plot(points_x,var.loc['Unconditioned'].values, linewidth = size_points*2, c = 'gray')
# axes[0,1].set_ylim(1E-18,1E-13)
axes[0,1].set_ylabel('Conditional Variance', fontsize = size_label)
axes[0,1].set_yscale('log')
axes[0,1].tick_params(axis = "y", which = 'both', right = True, direction = "in", labelsize = size_label)
axes[0,1].xaxis.set_minor_locator(AutoMinorLocator())
axes[0,1].tick_params(axis = "x", which = 'both', top = True, direction = "in", labelsize = size_label)

### SKEWNESS
# axes[1,0].text(0.9, 0.85, '(c)', transform=axes[1,0].transAxes, fontsize = size_label+constant, weight='bold')
# axes[1,0].scatter(points_x,skew.loc['Radius'].values, s = size_points, c = 'black')
# axes[1,0].scatter(points_x,skew.loc['Porosity'].values, s = size_points, c = 'purple')
# axes[1,0].scatter(points_x,skew.loc['Pressure'].values, s = size_points, c = 'red')
# axes[1,0].scatter(points_x,skew.loc['Tortuosity'].values, s = size_points, c = 'blue')
# axes[1,0].scatter(points_x,skew.loc['T'].values, s = size_points, c = 'green')
# axes[1,0].plot(points_x,skew.loc['Unconditioned'].values, linewidth = size_points*2, c = 'gray')
# axes[1,0].set_xlim(0,1)
# axes[1,0].set_ylim(0,8)
# axes[1,0].set_xlabel('Normalized Parameter', fontsize = size_label)
# axes[1,0].set_ylabel('Conditional Skewness', fontsize = size_label)
# axes[1,0].tick_params(axis = "y", which = 'both', right = True, direction = "in", labelsize = size_label)
# axes[1,0].xaxis.set_minor_locator(AutoMinorLocator())
# axes[1,0].tick_params(axis = "x", which = 'both', top = True, direction = "in", labelsize = size_label)

# ### KURTOSIS
# axes[1,1].text(0.9, 0.85, '(d)', transform=axes[1,1].transAxes, fontsize = size_label+constant, weight='bold')
# axes[1,1].scatter(points_x,kurt.loc['Radius'].values, s = size_points, label = r'$r_0$', c = 'black')
# axes[1,1].scatter(points_x,kurt.loc['Porosity'].values, s = size_points, label = r'$\phi_0$', c = 'purple')
# axes[1,1].scatter(points_x,kurt.loc['Pressure'].values, s = size_points, label = r'$p$', c = 'red')
# axes[1,1].scatter(points_x,kurt.loc['Tortuosity'].values, s = size_points, label = r'$\tau$', c = 'blue')
# axes[1,1].scatter(points_x,kurt.loc['T'].values, s = size_points, label = r'$T$', c = 'green')
# axes[1,1].plot(points_x,kurt.loc['Unconditioned'].values, linewidth = size_points*2, c = 'gray', label = r'$SM_Y$')
# axes[1,1].set_xlim(0,1)
# axes[1,1].set_ylim(1,100)
# axes[1,1].set_xlabel('Normalized Parameter', fontsize = size_label)
# axes[1,1].set_ylabel('Conditional Kurtosis', fontsize = size_label)
# axes[1,1].set_yscale('log')
# axes[1,1].tick_params(axis = "y", which = 'both', right = True, direction = "in", labelsize = size_label)
# axes[1,1].xaxis.set_minor_locator(AutoMinorLocator())
# axes[1,1].tick_params(axis = "x", which = 'both', top = True, direction = "in", labelsize = size_label)

plt.legend(bbox_to_anchor=(0.7, -0.15), scatterpoints = 3, frameon = False, prop={'size': 8}, ncol = 6)
plt.show()
# plt.savefig('conditional_mean.svg', format = 'svg')