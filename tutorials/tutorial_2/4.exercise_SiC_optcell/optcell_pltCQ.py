#!/bin/python

import matplotlib.pylab as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter

dft_energy = np.loadtxt('tmp.energy')
dft_rforce = np.loadtxt('tmp.stress')

fig, ax1 = plt.subplots(figsize=(8,6))

color1 = 'tab:blue'
ax1.plot(dft_energy,"o-", color=color1, label="Energy")
ax1.set_xlabel('# iteration',size=12)
ax1.set_ylabel('DFT Energy [Ha]', color=color1, size=12)
ax1.tick_params(axis='y', labelcolor=color1)
#ax1.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
ax1.grid(True)

color2 = 'tab:red'
ax2 = ax1.twinx()
ax2.plot(dft_rforce, 's--', color=color2, label='Stress')
ax2.set_ylabel('Stress Residual [GPa]', color=color2, size=12)
ax2.tick_params(axis='y', labelcolor=color2)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

fig.savefig('tmp.png', dpi=300,transparent=True)
plt.show()
