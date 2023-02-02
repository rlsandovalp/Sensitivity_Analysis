import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Information",
    page_icon="🔢",
    layout="wide",
    initial_sidebar_state="expanded"
)

r'''
## Information required to perform numerically a GSA

In order to perform a Global Sensitivity Analysis one needs to perform model evaluations of several 
combinations of the model parameters. The values of such parameters are typically selected randomly from 
their probability density function (pdf). 

In cases in which the information of the parameters of the analyzed model is limited, one can define upper and lower 
boundaries for the values of a parameter and then consider that the pdf of the model parameters is uniform 
(i.e., one considers that all parameter values within the identified range of variability are equally probable). 
However, in cases in which the pdf of the model parameters is known, one can sample ramdom values of the parameters from 
such a pdf.

At the end of the model evaluation process one should have a table in which for each model parameter combination one has 
the value of the model output. For example, if the model is 

$$
y = x_1^2 + 3x_2
$$

and the ranges of variability of $x_1$ and $x_2$ are (2,3) and (1,5), respectively; the model evaluation table could be:
'''

x1 = np.random.rand(10)+2
x2 = np.random.rand(10)*4+1
y = x1**2+3*x2

data = np.stack((x1,x2,y), axis = 1)

df = pd.DataFrame(data, columns=[r'x_1', r'x_2', r'y'])

st.table(df)

'''
Recall that the value of such combinations is random, thus, it is ok if you have different values. Note that the number of 
model evaluations in this table is only 10. Note that typically one requires many more to accurately perform Global Sensitivity 
Analysis. The number of required model evaluations grows (not linearly) with the number of model parameters and with the 
complexity of the model surface response. If the evaluation of the analyzed model is computationally expensive, one can perform 
the sensitivity analysis based on the outputs of a surrogate model.
'''