import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="Conditional Plots",
    page_icon="🆕",
    layout="wide",
    initial_sidebar_state="expanded"
)

def cond_plots(cond_x, cond_y, uncond):
    fig, ax = plt.subplots(1,2, figsize = (8, 3), sharey=True)
    for var in range(N_vars):
        if SM == 'Mean':
            UC = uncond[0]
            ax[var].plot(cond_x[var], cond_y[var], '.', label = '$E[y|x_'+str(var+1)+']$')
        else:
            UC = uncond[1]
            ax[var].plot(cond_x[var], cond_y[var], '.', label = '$V[y|x_'+str(var+1)+']$')
        ax[var].plot(np.linspace(xlims[var][0], xlims[var][1], N_intervals), UC*np.ones(N_intervals), label = 'Uncond')
        ax[var].set_xlim(xlims[var][0], xlims[var][1])
        ax[var].set_xlabel(vars[var])
        ax[var].legend()
    ax[0].set_ylabel(r'$y$')
    return fig

r'''
## Conditional Plots

The information collected to perform GSA can also be employed to construct conditional plots, i.e., plots that allow 
understanding the relationships between values of model parameters and statistical moments of model outputs. To construct 
this type of plot for the model parameter $x_i$ and the statistical moment $M$, one should evaluate $M(y)$ conditional to 
a value $x_i^v$ belonging to the support ($\Gamma_i$) of the model parameter (i.e., belonging to $\Gamma_i = [x_{i,min},x_{i,max}]$). 
For this one should evaluate

$$
M[y|x_i = x_i^v].
$$

In practice, since the number of model evaluations is limited because one can not evaluate the model an infinite number of times 
one evaluates the statistical moment of $y$ in the interval [$x_i^v-\varepsilon, x_i^v+\varepsilon$]. Note that smaller 
values of $\varepsilon$ are associated with more accurate representations of the conditional plots.

Thus, for our previous example $(y = x_1^2 + 3x_2)$ the conditional mean plot of $x_1$ is:
'''


N_sims = st.sidebar.selectbox("Number of evaluatiations of the model", [1000, 10000, 100000, 1000000], index = 0)
N_intervals = st.sidebar.selectbox("Number of intervals", [10, 20, 50, 100], index = 0)
SM = st.sidebar.selectbox('Statistical Moment', ['Mean', 'Variance'], index=0)


vars = ['x1', 'x2']
N_vars = len(vars)

cond_x = [[] for _ in range(N_vars)]
cond_y = [[] for _ in range(N_vars)]

xlims = [[2,3],[1,5]]
dx = [(xlims[i][1]-xlims[i][0])/N_intervals for i in range(N_vars)]

x = [np.random.rand(N_sims)*(xlims[i][1]-xlims[i][0])+xlims[i][0] for i in range(N_vars)]
y = x[0]**2+3*x[1]

uncond = np.mean(y)

uncond = [np.mean(y), np.var(y)]


for i in range(N_intervals):
    for var in range(N_vars):
        cond_x[var].append(xlims[var][0]+dx[var]*i+dx[var]/2)
        if SM == 'Mean':
            cond_y[var].append(np.mean(y[(x[var]>xlims[var][0]+dx[var]*i) & (x[var]<xlims[var][0]+dx[var]*(i+1))]))
        else:
            cond_y[var].append(np.var(y[(x[var]>xlims[var][0]+dx[var]*i) & (x[var]<xlims[var][0]+dx[var]*(i+1))]))



st.pyplot(cond_plots(cond_x, cond_y, uncond))