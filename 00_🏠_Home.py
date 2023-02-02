import streamlit as st

st.set_page_config(
    page_title="Sensitivity Analysis",
    page_icon="📈",
    layout="wide",
    initial_sidebar_state="expanded"
)

r"""
## Global Sensitivity Analysis

## Objective
To illustrate:
- What Sensitivity Analysis is and what Moment-based AMA indices are
- What information one needs to perform a numerical Sensitivity Analysis via AMA indices
- What additional products you can get with the information required for SA
- Example 1 - Ishigami Function
- Example 2 - Gas flow in low permeable materials (Under construction)

The theory and examples are taken from:
- Dell'Oca, A., Riva, M. and Guadagnini, A., 2017. Moment-based metrics for global sensitivity analysis of hydrological systems. Hydrology and Earth System Sciences, 21(12), pp.6219-6234. https://doi.org/10.5194/hess-21-6219-2017
- Sandoval, L., Riva, M., Colombo, I. and Guadagnini, A., 2022. Sensitivity Analysis and Quantification of the Role of Governing Transport Mechanisms and Parameters in a Gas Flow Model for Low-Permeability Porous Media. Transport in Porous Media, 142(3), pp.509-530. https://doi.org/10.1007/s11242-022-01755-x
"""