#import matplotlib.pyplot as plt
#import numpy as np
#from matplotlib.ticker import FuncFormatter

import numpy as np
import matplotlib
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


x = np.arange(4)
totals = [0.389, 0.274, 0.193, 0.143]
mkl = [0.389, 0.226, 0.0, 0.0]

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
ax.yaxis.set_major_formatter(formatter)
ax.set_ylabel('Contribution to total runtime')
ax.set_xlabel('Function group')
ax.set_yticks(np.arange(0.1, 0.5, step=0.1))
plt.bar(x, totals, color='#e41a1c')
plt.bar(x, mkl, color='#377eb8')
plt.xticks(x, ('BLAS', 'ETC', 'USER', 'MPI'))
plt.legend(['Total', 'MKL'])
fig.set_size_inches(6, 4)
filename = 'MnGe_6x6x6_crayreport.pdf' 
fig.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')
#plt.show()
