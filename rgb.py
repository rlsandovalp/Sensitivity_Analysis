import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

######################## Input parameters #########################
files = 100
distribution = 'lognormal'       #uniform, lognormal, normal
nclass = 100

######################## Initialize processing files #########################
param_list = {'Radius': 0,'Porosity': 1,'Pressure': 2,'Tortuosity': 3,'q': 4,'t': 5,'Overburden Pressure': 6,'k': 7, 'Fractal Dimension': 8, 'Isosteric heat': 9, 'lp0': 10, 'a0': 11, 'a1': 12, 'beta': 13, 'T': 14}
Npar = 15
hf = h5py.File(distribution+"/results_wu_"+distribution+"1.hdf5", 'r')
Npoints = (np.array((hf.get('Coordinates')).get('Radius'))).size
hf.close()
Sampling_Points = np.zeros((Npoints*files,Npar))
Results = np.zeros((Npoints*files,3))

######################## Read files and join them in only one array #########################
for j in range (0,files):
    filename = distribution+"/results_wu_"+distribution+str(j)+".hdf5"
    hf = h5py.File(filename, 'r')
    for id_par,item_par in enumerate(param_list):
        Sampling_Points[j*Npoints:(j+1)*Npoints,id_par] = np.array((hf.get('Coordinates')).get(item_par))
    Results[j*Npoints:(j+1)*Npoints,:] = np.array((hf.get('Flow')).get('Results'))[:,(0,1,2)]
    hf.close()

######################## Mean #########################

Processing = np.zeros((Npoints*files,5))
Processing[:,2:5] = Results[:,0:3]
Output = np.zeros((nclass,nclass,3))

fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (10.5,9))

print('Figure 1')
Processing[:,0:2] = Sampling_Points[:,(param_list['Pressure'],param_list['Radius'])]
x_min = Processing[:,0].min()
x_max = Processing[:,0].max()
y_min = Processing[:,1].min()
y_max = Processing[:,1].max()
discret_x = (x_max-x_min)/nclass
discret_y = (y_max-y_min)/nclass

for iclass1 in range (nclass):
    processing_class = Processing[Processing[:,1]<y_min+discret_y*(iclass1+1)]
    processing_class = processing_class[processing_class[:,1]>y_min+discret_y*iclass1]
    for iclass2 in range (0,nclass):
        processing_class1 = processing_class[processing_class[:,0]<x_min+discret_x*(iclass2+1)]
        processing_class1 = processing_class1[processing_class1[:,0]>x_min+discret_x*iclass2]
        Output[iclass1,iclass2,:] = np.mean(processing_class1[:,2:5], axis = 0)

axes[0,0].imshow(Output, origin = 'lower', extent = [0,100,0,100])
axes[0,0].set_xticks(np.arange(0, 101, 50))
axes[0,0].set_yticks(np.arange(0, 101, 50))
axes[0,0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[0,0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[0,0].tick_params(axis = "x", labelsize = 16, top = 1)
axes[0,0].tick_params(axis = "y", labelsize = 16)
axes[0,0].xaxis.set_ticks_position('top')
axes[0,0].xaxis.set_label_position('top')
axes[0,0].set_xlabel(r'$p$', fontsize = 18)
axes[0,0].set_ylabel(r'$r_o$', fontsize = 18)


print('Figure 2')
Processing[:,0:2] = Sampling_Points[:,(param_list['Pressure'],param_list['Porosity'])]
x_min = Processing[:,0].min()
x_max = Processing[:,0].max()
y_min = Processing[:,1].min()
y_max = Processing[:,1].max()
discret_x = (x_max-x_min)/nclass
discret_y = (y_max-y_min)/nclass

for iclass1 in range (nclass):
    processing_class = Processing[Processing[:,1]<y_min+discret_y*(iclass1+1)]
    processing_class = processing_class[processing_class[:,1]>y_min+discret_y*iclass1]
    for iclass2 in range (0,nclass):
        processing_class1 = processing_class[processing_class[:,0]<x_min+discret_x*(iclass2+1)]
        processing_class1 = processing_class1[processing_class1[:,0]>x_min+discret_x*iclass2]
        Output[iclass1,iclass2,:] = np.mean(processing_class1[:,2:5], axis = 0)

axes[0,1].imshow(Output, origin = 'lower', extent = [0,100,0,100])
axes[0,1].set_xticks(np.arange(0, 101, 50))
axes[0,1].set_yticks(np.arange(0, 101, 50))
axes[0,1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[0,1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[0,1].tick_params(axis = "x", labelsize = 16)
axes[0,1].tick_params(axis = "y", labelsize = 16)
axes[0,1].xaxis.set_ticks_position('top')
axes[0,1].xaxis.set_label_position('top')
axes[0,1].set_xlabel(r'$p$', fontsize = 18)
axes[0,1].set_ylabel(r'$\phi_o$', fontsize = 18)



print('Figure 3')
Processing[:,0:2] = Sampling_Points[:,(param_list['Radius'],param_list['Porosity'])]
x_min = Processing[:,0].min()
x_max = Processing[:,0].max()
y_min = Processing[:,1].min()
y_max = Processing[:,1].max()
discret_x = (x_max-x_min)/nclass
discret_y = (y_max-y_min)/nclass

for iclass1 in range (nclass):
    processing_class = Processing[Processing[:,1]<y_min+discret_y*(iclass1+1)]
    processing_class = processing_class[processing_class[:,1]>y_min+discret_y*iclass1]
    for iclass2 in range (0,nclass):
        processing_class1 = processing_class[processing_class[:,0]<x_min+discret_x*(iclass2+1)]
        processing_class1 = processing_class1[processing_class1[:,0]>x_min+discret_x*iclass2]
        Output[iclass1,iclass2,:] = np.mean(processing_class1[:,2:5], axis = 0)

axes[1,0].imshow(Output, origin = 'lower', extent = [0,100,0,100])
axes[1,0].set_xticks(np.arange(0, 101, 50))
axes[1,0].set_yticks(np.arange(0, 101, 50))
axes[1,0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[1,0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[1,0].tick_params(axis = "x", labelsize = 16)
axes[1,0].tick_params(axis = "y", labelsize = 16)
axes[1,0].set_xlabel(r'$r_o$', fontsize = 18)
axes[1,0].set_ylabel(r'$\phi_o$', fontsize = 18)



print('Figure 4')
Processing[:,0:2] = Sampling_Points[:,(param_list['Pressure'],param_list['T'])]
x_min = Processing[:,0].min()
x_max = Processing[:,0].max()
y_min = Processing[:,1].min()
y_max = Processing[:,1].max()
discret_x = (x_max-x_min)/nclass
discret_y = (y_max-y_min)/nclass

for iclass1 in range (nclass):
    processing_class = Processing[Processing[:,1]<y_min+discret_y*(iclass1+1)]
    processing_class = processing_class[processing_class[:,1]>y_min+discret_y*iclass1]
    for iclass2 in range (0,nclass):
        processing_class1 = processing_class[processing_class[:,0]<x_min+discret_x*(iclass2+1)]
        processing_class1 = processing_class1[processing_class1[:,0]>x_min+discret_x*iclass2]
        Output[iclass1,iclass2,:] = np.mean(processing_class1[:,2:5], axis = 0)

axes[1,1].imshow(Output, origin = 'lower', extent = [0,100,0,100])
axes[1,1].set_xticks(np.arange(0, 101, 50))
axes[1,1].set_yticks(np.arange(0, 101, 50))
axes[1,1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[1,1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
axes[1,1].tick_params(axis = "x", labelsize = 16)
axes[1,1].tick_params(axis = "y", labelsize = 16)
axes[1,1].set_xlabel(r'$p$', fontsize = 18)
axes[1,1].set_ylabel(r'$T$', fontsize = 18)
plt.tight_layout()
plt.show()