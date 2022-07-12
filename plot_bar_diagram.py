import numpy as np
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('../SA_results.txt', sep = '\t', index_col = 0)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax3 = fig.add_axes([0.5, 0.3, 0.35, 0.5])
data.AMAVi.plot(kind = 'bar', color = 'tab:blue', ax = ax, width = 0.4, position = 0)
data.Si.plot(kind = 'bar', color = 'tab:red', ax = ax2, width = 0.4, position = 1)
data.iloc[0:6].AMAEi.plot(kind = 'bar', color='tab:green', ax = ax3)
ax.set_ylabel(r'AMA$V_{x_i}$')
ax2.set_ylabel(r'$S_{x_i}$')
ax3.set_ylabel(r'AMA$E_{x_i}$')

plt.show()