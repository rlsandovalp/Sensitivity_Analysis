import numpy as np
import pandas as pd
import math
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(
    page_title="Ishigami Function",
    page_icon="⛩️",
    layout="wide",
    initial_sidebar_state="expanded"
)

def cond_plots(cond_x, cond_y, uncond):
    fig, ax = plt.subplots(1,3, figsize = (8, 3), sharey=True)
    for var in range(N_vars):
        if SM == 'Mean':
            UC = uncond[0]
            ax[var].plot(cond_x[var], cond_y[var], '.', label = '$E[y|x_'+str(var+1)+']$')
            ax[var].plot(cond_x[var], cond_y_ana_mean[var], label = 'Analytical')
        else:
            UC = uncond[1]
            ax[var].plot(cond_x[var], cond_y[var], '.', label = '$V[y|x_'+str(var+1)+']$')
            ax[var].plot(cond_x[var], cond_y_ana_var[var], label = 'Analytical')
        ax[var].plot(np.linspace(xlims[var][0], xlims[var][1], N_intervals), UC*np.ones(N_intervals), label = 'Uncond')
        ax[var].set_xlim(xlims[var][0], xlims[var][1])
        ax[var].set_xlabel(vars[var])
        ax[var].legend()
    ax[0].set_ylabel(r'$y$')
    return fig

def AMAE_indices(AMA):
    x = ('x1', 'x2', 'x3')
    if SM == 'Mean':
        aronne = [0.75, 0.64, 0.0]
    else:
        aronne = [0.40, 0.29, 0.84]
    X_axis = np.arange(len(x))
    fig, ax = plt.subplots(figsize = (8,3))
    ax.bar(X_axis - 0.2, AMA, width = 0.4, label = 'Numerical')
    ax.bar(X_axis + 0.2, aronne, width = 0.4, label = 'Analytical')
    ax.set_xticks(X_axis, x)
    ax.legend()

    return fig
r'''
## Ishigami Function

The nonlinear and non-monotonic Ishigami function

$$
y(\mathbf{x}) = ISH(\mathbf{x}) = \sin(2\pi x_1-\pi)+a\sin^2(2\pi x_2-pi)+b(2\pi x_3-\pi)^4\sin(2\pi x_1-\pi)
$$

is widely used in the literature (e.g., Homma and Saltelli, 1996; Chun et al., 2000; Borgonovo, 2007; Sudret, 2008; 
Crestaux et al., 2009; Borgonovo et al. 2011) to benchmark GSA methods. Here, $x_i (i = 1, 2, 3)$ are i.i.d. random 
variables uniformly distributed within the interval [0, 1]. Unconditional mean E[ISH], variance, V[ISH], skewness, γ[ISH],
and kurtosis, $k$[ISH], can be evaluated analytically as reported in **Dell'Oca et al. (2017). The study from which this 
example is taken**.

Here we are going to approximate numerically to the evaluation of the moment-based AMA indices and to the 
variance-based Sobol indices. First, we can construct the conditional plots.

### Conditional Plots
'''

a, b = 5, 0.1

N_sims = st.sidebar.selectbox("Number of evaluatiations of the model", [1000, 10000, 100000, 1000000], index = 0)
N_intervals = st.sidebar.selectbox("Number of intervals", [10, 20, 50, 100], index = 0)
SM = st.sidebar.selectbox('Statistical Moment', ['Mean', 'Variance'], index=0)

vars = ['x1', 'x2', 'x3']
N_vars = len(vars)

cond_x = [[] for _ in range(N_vars)]
cond_y = [[] for _ in range(N_vars)]

xlims = [[0,1],[0,1],[0,1]]

dx = [(xlims[i][1]-xlims[i][0])/N_intervals for i in range(N_vars)]

x = [np.random.rand(N_sims) for _ in range(N_vars)]
y = np.sin(2*math.pi*x[0]-math.pi) + a*np.sin(2*math.pi*x[1]-math.pi)**2 + b*(2*math.pi*x[2]-math.pi)**4*np.sin(2*math.pi*x[0]-math.pi)

cond_y_ana_mean = []
cond_y_ana_var = []

uncond = [np.mean(y), np.var(y)]

for i in range(N_intervals):
    for var in range(N_vars):
        cond_x[var].append(xlims[var][0]+dx[var]*i+dx[var]/2)
        if SM == 'Mean':
            cond_y[var].append(np.mean(y[(x[var]>xlims[var][0]+dx[var]*i) & (x[var]<xlims[var][0]+dx[var]*(i+1))]))
        else:
            cond_y[var].append(np.var(y[(x[var]>xlims[var][0]+dx[var]*i) & (x[var]<xlims[var][0]+dx[var]*(i+1))]))

cond_y_ana_mean.append(a/2 - 1/5*(5+b*math.pi**4)*np.sin(2*math.pi*np.array(cond_x[0])))
cond_y_ana_mean.append(a*np.sin(2*math.pi*np.array(cond_x[1]))**2)
cond_y_ana_mean.append(a/2*np.ones(N_intervals))

cond_y_ana_var.append(a**2/8 + 8*b**2*math.pi**8/225*(1-np.cos(4*math.pi*np.array(cond_x[0]))))
cond_y_ana_var.append(0.5+b*math.pi**4*(0.2+b/18*math.pi**4)*np.ones(N_intervals))
cond_y_ana_var.append(a**2/8+0.5*(1+b*math.pi**4*(1-2*np.array(cond_x[2]))**4)**2)


st.pyplot(cond_plots(cond_x, cond_y, uncond))

'''
Now, we can evaluate the AMA indices and compare them with the results reported in **Dell'Oca et al. (2017)**

'''
if SM == 'Mean':
    AMA = [np.mean(abs(uncond[0] - np.array(cond_y[i])))/uncond[0] for i in range(N_vars)]
else:
    AMA = [np.mean(abs(uncond[1] - np.array(cond_y[i])))/uncond[1] for i in range(N_vars)]

st.pyplot(AMAE_indices(AMA))