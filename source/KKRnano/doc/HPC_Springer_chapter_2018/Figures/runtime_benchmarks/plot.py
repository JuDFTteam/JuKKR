#import matplotlib.pyplot as plt
#import numpy as np
#from matplotlib.ticker import FuncFormatter

import numpy as np
import matplotlib
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

width = 0.2       # the width of the bars
x = np.arange(2)
total = [293, 371]
single_cell = [110, 107]
multi_scattering = [125, 122]
electrostatic = [3, 10]

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'

formatter = FuncFormatter(to_percent)

fig, ax = plt.subplots()
#ax.yaxis.set_major_formatter(formatter)
ax.set_xlabel('# Atoms')
ax.set_ylabel('Runtime (s)')
ax.set_yticks(np.arange(0, 400, step=50))
plt.bar(x, total, width, color='#e41a1c')
plt.bar(x+width, multi_scattering, width, color='#4daf4a')
plt.bar(x+2*width, single_cell, width, color='#377eb8')
plt.bar(x+3*width, electrostatic, width, color='#984ea3')
plt.xticks(x, ('1728', '13824'))
plt.legend(['Total', 'Single-cell', 'Multi-scattering', 'Electrostatics'])
fig.set_size_inches(7, 4)
filename = 'MnGe_benchmarks.pdf' 
fig.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')
#plt.show()
