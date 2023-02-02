import streamlit as st

st.set_page_config(
    page_title="Sensitivity Analysis and Moment-Based AMA indices",
    page_icon="💭",
    layout="wide",
    initial_sidebar_state="expanded"
)

r'''
## Global Sensitivity Analysis

Global Sensitivity Analysis (GSA) techniques are important tools enabling us to:
- Quantify uncertainty, 
- Enhance our understanding of the relationships between model inputs and outputs, and 
- Tackle the challenges of model- and data-driven design of experiments (Dell’Oca et al. 2017). 

Hence, GSA techniques may be effectively used in the context of parametrized models to:
- Quantify and rank the contribution of our lack of knowledge on each model parameter to the uncertainty associated with model outputs,
- Identify model input–output relationships, and 
- Enhance the quality of parameter estimation workflows, upon focusing efforts on parameters with the highest influence to target model outputs (Saltelli et al. 2010). 

In cases where parameters associated with a model have already been estimated (e.g., through model calibration), 
the main purpose of a GSA is to assist quantification of the uncertainty still 
remaining after model calibration, thus guiding additional efforts for its characterization 
(e.g., Dell’Oca et al. (2020) and references therein). 

The probability density function (pdf) related to each model parameter at this step might 
differ from the one employed before model calibration and some model parameters might be associated 
with a reduced uncertainty. In cases where processes are described through black-box 
models, GSA can be employed to quantify the influence that the variability 
of hyperparameters embedded in these models can have on their outcomes. We note that if uncertainty of 
some model parameters is further constrained, for example through (stochastic) inverse modeling 
(e.g., Ceresa et al. (2021)), results of the uncertainty quantification might also change. 

In this site we illustrate the methodological framework and the workflow required for GSA of parametrized 
models, provide the elements to perform such an analysis for diverse scenarios, and provide two applications 
GSA techniques.

## Global Sensitivity Analysis Approaches

### Variance-based Sobol’ Indices

Sobol' indices (Saltelli and Sobol’ 1995) can assist the appraisal and quantification of the 
relative expected reduction of the variance of a target model output due to knowledge 
of (or conditioning on) a given model parameter, which would otherwise be subject to 
uncertainty. In this context, considering a model output $y$, which depends on $N$ 
random parameters collected in vector $\textbf{x} = (x_1, x_2, ..., x_N)$ and 
defined within the space $\Gamma = \Gamma_1 \times \Gamma_2 \times ... \times \Gamma_N$ ($\Gamma_i = [x_{i,min},x_{i,max}]$ 
corresponding to the support of the $i$-th parameter, $x_i$), the principal Sobol' 
index $S_{x_i}$ associated with a given model parameter $x_i$ is evaluated as

$$
S_{x_i} = \frac{V\left[E\left[y|x_i\right]\right]}{V\left[y\right]}.
$$

Here, $E\left[\cdot\right]$ and $V\left[\cdot\right]$ represent expectation and 
variance operators, respectively; the notation $y|x_i$ denotes conditioning of $y$ 
on $x_i$. Note that $S_{x_i}$ describes the relative contribution to $V[y]$ due to 
variability of only $x_i$. Joint contributions of $x_i$ with other model parameters 
included in $\textbf{x}$ to the variance of $y$ are embedded in the total Sobol' 
indices (details not shown). We recall that relying on Sobol' indices to diagnose 
the relative importance of uncertain model parameters to model outputs is tantamount 
to identifying uncertainty with the concept of variance of a pdf. As such, while Sobol’ 
indices are characterized by a conceptual simplicity and straightforward implementation 
and use, they provide only limited information about the way variations of model 
parameters can influence the complete pdf of model outputs.

### Moment-Based AMA Indices

The recent moment-based GSA approach proposed by Dell’Oca et al. (2017, 2020) rests on the idea 
that the quantification of the effects of model parameter uncertainty on various statistical 
moments of the ensuing pdf of model outputs can provide enhanced understanding of model functioning. 
Dell’Oca et al. (2017) introduce Moment-Based sensitivity metrics (termed AMA indices after the 
initials of the authors) according to which one can evaluate the influence of uncertain model 
parameter on key elements of the model output pdf, as embedded in its associated statistical 
moments. The AMA indices are defined as follows (Dell’Oca et al. (2017)):

$$
\text{AMA}M_{x_i} = \frac{1}{\left|M[y]\right|} E\left[\left|M[y]-M\left[y|x_i\right]\right|\right]
$$

Here, AMA$M_{x_i}$ represents the indices associated with a model parameter ${x_i}$ and a 
given statistical moment $M$ of the pdf of model output $y$ (considering the first four 
statistical moments of $y$, $M=E$ for the mean, $M=V$ for the variance, $M=\gamma$ for the 
skewness, and $M=k$ for the kurtosis). The AMA indices are intended to quantify the expected 
change of each statistical moment of $y$ due to our knowledge of $x_i$. Large values of these 
indices indicate that variations of the associated parameter strongly affect the statistical 
moments of $y$.


## References

- Ceresa, L., Guadagnini, A., Porta, G.M., Riva, M.: Formulation and probabilistic assessment of reversible biodegradation pathway of Diclofenac in groundwater. Water Res. 204(117), 466 (2021). https://doi.org/10.1016/J.WATRES.2021.117466
- Dell’Oca, A., Riva, M., Guadagnini, A.: Moment-based metrics for global sensitivity analysis of hydrological systems. Hydrol. Earth Syst. Sci. 21(12), 6219–6234 (2017). https://doi.org/10.5194/hess-21-6219-2017
- Dell’Oca, A., Riva, M., Guadagnini, A.: Global sensitivity analysis for multiple interpretive models with uncertain parameters. Water Resour. Res. 56(2), 1–20 (2020). https://doi.org/10.1029/2019WR025754
- Saltelli, A., Sobol’, I.M.: Sensitivity analysis for nonlinear mathematical models: numerical experience (in Russian). Math. Models Comput. Exp. 7(11), 16–28 (1995)
- Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., Tarantola, S.: Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Comput. Phys. Commun. 181(2), 259–270 (2010)

'''


