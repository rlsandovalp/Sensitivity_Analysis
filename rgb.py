import numpy as np
import pandas as pd
import h5py
import os
import matplotlib.pyplot as plt

######################## Input parameters #########################
files = 2
distribution = 'lognormal'       #uniform, lognormal, normal
nclass = 100
# The following list must be filled with the pairs of variables for each figure
# figures = [[2,0],[1,0],[2,14]]
figures = [[2,0],[2,1],[2,3],[2,14],[1,0],[3,0]]
######################## Initialize processing files #########################
param_list = ['Radius','Porosity','Pressure','Tortuosity','q','t','Overburden Pressure','k', 'Fractal Dimension', 'Isosteric heat', 'lp0', 'a0', 'a1', 'beta', 'T']
Npar = np.size(param_list)
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
    print('I joined '+str(j+1)+' files')

######################## Mean #########################

for figura in figures:
    param_x = figura[0]
    param_y = figura[1]
    Processing = np.zeros((Npoints*files,5))
    Output = np.zeros((3,nclass,nclass))

    Processing[:,0:2] = Sampling_Points[:,(param_y,param_x)]
    Processing[:,2:5] = Results[:,0:3]

    x_min = Sampling_Points[:,param_x].min()
    x_max = Sampling_Points[:,param_x].max()
    y_min = Sampling_Points[:,param_y].min()
    y_max = Sampling_Points[:,param_y].max()

    discret_1 = (y_max-y_min)/nclass
    discret_2 = (x_max-x_min)/nclass

    for iclass1 in range (0,nclass):
        print(iclass1)
        processing_class = Processing[Processing[:,0]<y_min+discret_1*(iclass1+1)]
        processing_class = processing_class[processing_class[:,0]>y_min+discret_1*iclass1]
        for iclass2 in range (0,nclass):
            processing_class1 = processing_class[processing_class[:,1]<x_min+discret_2*(iclass2+1)]
            processing_class1 = processing_class1[processing_class1[:,1]>x_min+discret_2*iclass2]
            for rgb in range(0,3):
                Output[rgb,iclass1,iclass2] = np.mean(processing_class1[:,rgb+2])

    Hola = np.zeros((nclass,nclass,3))
    Hola[:,:,0] = Output[0,:,:]
    Hola[:,:,1] = Output[1,:,:]
    Hola[:,:,2] = Output[2,:,:]

    plt.imshow(Hola, origin = 'lower', extent = [0,1,0,1])
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.xlabel(param_list[param_x], fontsize = 14)
    plt.ylabel(param_list[param_y], fontsize = 14)
    # plt.savefig('../../Figures/'+str(param_list[param_x])+'_'+str(param_list[param_y])+'_'+distribution+'.png', dpi=300)


