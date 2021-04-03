import numpy as np
import h5py
from scipy.stats import skew, norm, kurtosis
from matplotlib import pyplot as plt 
import matplotlib.ticker as mtick
from timeit import default_timer as timer

######################## Input parameters #########################
files = 1000
distribution = 'lognormal'            #uniform, lognormal, normal
old = 0               # 1 Yes, 2 No
nclass = 100


######################## Initialize processing files #########################
param_list = ['Radius','Porosity','Pressure','Tortuosity','q','t','Overburden Pressure','k', 'Fractal Dimension', 'Isosteric heat', 'lp0', 'a0', 'a1', 'beta', 'T']
Npar = np.size(param_list)
hf = h5py.File(distribution+"/results_wu_"+distribution+"1.hdf5", 'r')
if old == 1:
    hf.close()
    hf = h5py.File(distribution+"_old/results_wu_"+distribution+"1.hdf5", 'r')

Npoints = (np.array((hf.get('Coordinates')).get('Radius'))).size
hf.close()
Sampling_Points = np.zeros((Npoints*files,Npar))
Results = np.zeros(Npoints*files)

######################## Read files and join them in only one array #########################
for j in range (0,files):
    if old == 1:
        filename = distribution+"_old/results_wu_"+distribution+str(j)+".hdf5"
    else:
        filename = distribution+"/results_wu_"+distribution+str(j)+".hdf5"
    hf = h5py.File(filename, 'r')
    for id_par,item_par in enumerate(param_list):
        Sampling_Points[j*Npoints:(j+1)*Npoints,id_par] = np.array((hf.get('Coordinates')).get(item_par))
    if old == 1:
        Results[j*Npoints:(j+1)*Npoints] = np.array((hf.get('Flow')).get('J_0.1'))[:,0]
    else:
        Results[j*Npoints:(j+1)*Npoints] = np.array((hf.get('Flow')).get('Results'))[:,(3)]
    hf.close()
    print(j)

print('Mean ', np.mean(Results))
print('Var ', np.var(Results))
print('Skew ', skew(Results))
print('Kurt ', 3 + kurtosis(Results))

cond_mean=np.zeros((Npar,nclass))
cond_var=np.zeros((Npar,nclass))
cond_skew=np.zeros((Npar,nclass))
cond_kurt=np.zeros((Npar,nclass))

AMAE=np.zeros(Npar)
AMAV=np.zeros(Npar)
Si=np.zeros(Npar)
AMAskew=np.zeros(Npar)
AMAkurt=np.zeros(Npar)

Output = np.zeros((Npoints*files,2))

for ipar,item_par in enumerate (param_list):
    print(item_par)
    Output[:,0] = Sampling_Points[:,ipar]
    par_min = Sampling_Points[:,ipar].min()
    par_max = Sampling_Points[:,ipar].max()
    Output[:,1] = Results
    discret = (par_max-par_min)/nclass
    for iclass in range (0,nclass):
        output_class = (Output[(Output[:,0]<par_min+discret*(iclass+1)) & (Output[:,0]>par_min+discret*iclass)])[:,1]
        cond_mean[ipar,iclass]=np.mean(output_class)
        cond_var[ipar,iclass]=np.var(output_class)
        cond_skew[ipar,iclass]=skew(np.reshape(output_class,-1))
        cond_kurt[ipar,iclass]=3+kurtosis(np.reshape(output_class,-1))
    Si[ipar]=(np.var(cond_mean[ipar,:],ddof=1))/np.var(Output[:,1])
    AMAE[ipar]=np.mean(abs(cond_mean[ipar,:]-np.mean(Output[:,1])))/np.mean(Output[:,1])
    AMAV[ipar]=np.mean(abs(cond_var[ipar,:]-np.var(Output[:,1])))/np.var(Output[:,1])
    AMAskew[ipar]=np.mean(abs(cond_skew[ipar,:]-skew(np.reshape(Output[:,1],-1))))/skew(np.reshape(Output[:,1],-1))
    AMAkurt[ipar]=np.mean(abs(cond_kurt[ipar,:]-(3+kurtosis(np.reshape(Output[:,1],-1)))))/(3+kurtosis(np.reshape(Output[:,1],-1)))

np.savetxt('SA_results/'+distribution+'_AMA.txt',(AMAE,AMAV,AMAskew,AMAkurt,Si))


# Comment/uncomment these lines if you want to save the conditioned values of your statistical moments.
np.savetxt('SA_results/'+distribution+'_CM.txt',cond_mean)
np.savetxt('SA_results/'+distribution+'_CV.txt',cond_var)
np.savetxt('SA_results/'+distribution+'_CS.txt',cond_skew)
np.savetxt('SA_results/'+distribution+'_CK.txt',cond_kurt)