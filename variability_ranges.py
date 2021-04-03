import numpy as np

points_per_interval = 1000

# Define ranges of variability of the parameters
p_range = np.linspace (500_000, 50_000_000, points_per_interval)          # Pore pressure                         Pa
ro_range = np.linspace(2E-9, 1E-7, points_per_interval)                # Pore radius                           m
por_range = np.linspace(0.005, 0.1, points_per_interval)               # Porosity                              -
tor_range = np.linspace(2.8, 5.8, points_per_interval)                 # Tortuosity                            -
op_range = np.linspace(51_000_000, 90_000_000, points_per_interval)        # Overburden pressure                   Pa
q_range = np.linspace(0.014, 0.056, points_per_interval)               # q exponent porosity power law         -
t_range = np.linspace(0.02, 0.04, points_per_interval)                 # t exponent pore radius power law      -
k_range = np.linspace(0.1, 2, points_per_interval)                     # Blockage/Migration ratio              -
fd_range = np.linspace(2.1, 2.9, points_per_interval)                  # Fractal dimension of pore wall        -
ih_range = np.linspace(12000, 16000, points_per_interval)              # Fractal dimension of pore wall        -
lp0_range = np.linspace(41_000_000, 128_000_000, points_per_interval)      # Langmuir pressure                     Pa
a0_range = np.linspace(1.02, 1.36, points_per_interval)                # Alpha zero                            -
a1_range = np.linspace(2, 6, points_per_interval)                      # Alpha one                             -
beta_range = np.linspace(0.2, 0.6, points_per_interval)                # Beta                                  -
T_range = np.linspace(337, 473, points_per_interval)                   # Temperature                           K  https://www.sciencedirect.com/science/article/pii/0378381295027716