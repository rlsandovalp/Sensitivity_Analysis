import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis

######################## Input file #########################
source = 'MC_S2_out.csv' 
nclass = 50

######################## Initialize processing files #########################
param_list = ['k_over','k_fault','time']
Npar = np.size(param_list)
table = pd.read_csv(source).values
Sampling_Points = table[:,0:Npar]
Results = table[:,Npar]

######################## Print unconditional statistical moments #########################
print('Mean ', np.mean(Results))
print('Var ', np.var(Results))
print('Skew ', skew(Results))
print('Kurt ', 3 + kurtosis(Results))

######################## Create matrices to save the conditional moments #########################
cond_mean=np.zeros((Npar,nclass))
cond_var=np.zeros((Npar,nclass))
cond_skew=np.zeros((Npar,nclass))
cond_kurt=np.zeros((Npar,nclass))

######################## Create matrices to save the GSA indices (AMA and sobol) #########################

AMAE=np.zeros(Npar)
AMAV=np.zeros(Npar)
Si=np.zeros(Npar)
AMAskew=np.zeros(Npar)
AMAkurt=np.zeros(Npar)

Output = np.zeros((np.size(Results),2))


######################## Do the GSA #########################
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

np.savetxt('results_SA.txt',(AMAE,AMAV,AMAskew,AMAkurt,Si))

# Comment/uncomment these lines if you want to save the conditioned values of your statistical moments.
np.savetxt('CM.txt',cond_mean)
np.savetxt('CV.txt',cond_var)
np.savetxt('CS.txt',cond_skew)
np.savetxt('CK.txt',cond_kurt)