# Link for skripts
# http://matplotlib.org/gallery.html

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import brewer2mpl

Ryd2mev = 13.605698066*1000.0
ev2mev = 1000.0

# activate LATEX font
plt.rc('text', usetex=True)

# Fonts
font = {'family': 'arial',
        'color':  'black',
        'weight': 'normal',
        'size': 9,
        }


# Load the DOS files - KKRnano
energy_socscale_ferro90             = np.loadtxt("energies_ferro_MnGe.txt",comments=['#','&'])
energy_socscale_bp                  = np.loadtxt("energies_bp_MnGe.txt",comments=['#','&'])
#energy_socscale_3q                  = np.loadtxt("energy_socscale_3q",comments=['#','&'])
moments_socscale_ferro90             = np.loadtxt("moments_ferro_MnGe.txt",comments=['#','&'])
moments_socscale_bp                  = np.loadtxt("moments_bp_MnGe.txt",comments=['#','&'])
#moments_socscale_3q                  = np.loadtxt("moments_socscale_3q",comments=['#','&'])

x       =  energy_socscale_ferro90[:,0]
energy_socscale_bp[:,1] = (energy_socscale_bp[:,1] - energy_socscale_ferro90[:,1])*ev2mev/864.0
#energy_socscale_3q[:,1] = (energy_socscale_3q[:,1] - energy_socscale_ferro90[:,1])*ev2mev/864.0
energy_socscale_ferro90[:,1] = (energy_socscale_ferro90[:,1] - energy_socscale_ferro90[:,1])*ev2mev/864.0

moments_socscale_ferro90[:,1] = (moments_socscale_ferro90[:,1])/864.0
moments_socscale_bp[:,1] = (moments_socscale_bp[:,1])/864.0
#moments_socscale_3q[:,1] = (moments_socscale_3q[:,1])/864.0



# add l=s,p,d,f contributions
plot_ingredients  = []
plot_ingredients.append(energy_socscale_ferro90[:,1])  
plot_ingredients.append(energy_socscale_bp[:,1])  
#plot_ingredients.append(energy_socscale_3q[:,1])  
plot_ingredients.append(moments_socscale_ferro90[:,1])  
plot_ingredients.append(moments_socscale_bp[:,1])  
#plot_ingredients.append(moments_socscale_3q[:,1])  



########################################################################
# Set colorbrewer colors for the plot
########################################################################
set2 = brewer2mpl.get_map('Set1', 'qualitative', 4).mpl_colors

########################################################################
# Plot files
########################################################################

# single plot

# multiple subplots
f, axarr = plt.subplots(2, sharex=True) 

pos  = ['Ferromagnet','Bloch point','3Q']

x_lim = [0.95,4.05]
y_lim = [-5,1]
y_lim2 = [2.05,2.20]

axarr[0].plot([1.0,5.0], [0,0], linestyle='--', color='black', linewidth=1.2)
#axarr[1].plot([4.8,4.8], [0.6,3.2], linestyle='--', color='black', linewidth=1.2)
#axarr[1].plot([4.5,5.2], [2.0,2.0], linestyle='--', color='black', linewidth=1.2)
#axarr[0].annotate('exp.', xy=(4.77, -5.0), xycoords='data', rotation=90)
#axarr[1].annotate('exp.', xy=(4.77, 1.3), xycoords='data', rotation=90)
#axarr[1].annotate('exp.', xy=(5.05, 1.75), xycoords='data')
#plt.plot([-10,10], [0,0], linestyle='-', color='black', linewidth=0.8)
#plt.plot([0,0], [-10,10], linestyle='-', color='black', linewidth=0.8)
#plt.plot([-10,10], [0,0], linestyle='-', color='black', linewidth=0.8)

#axarr[0].plot(energy_socscale_ferro90[:,0], energy_socscale_ferro90[:,1], marker='o', ls='-', c=set2[0],  lw=1.0, label=pos[0] )
axarr[0].plot(energy_socscale_bp[:,0], energy_socscale_bp[:,1], marker='o', ls='-', c=set2[1],  lw=1.0, label=pos[1] )
#axarr[0].plot(energy_socscale_3q[:,0], energy_socscale_3q[:,1], marker='o', ls='-', c=set2[2],  lw=1.0, label=pos[2] )
axarr[1].plot(moments_socscale_ferro90[:,0], moments_socscale_ferro90[:,1], marker='o', ls='-', c=set2[0],  lw=1.0, label=pos[0] )
axarr[1].plot(moments_socscale_bp[:,0], moments_socscale_bp[:,1], marker='o', ls='-', c=set2[1],  lw=1.0, label=pos[1] )
#axarr[1].plot(moments_socscale_3q[:,0], moments_socscale_3q[:,1], marker='o', ls='-', c=set2[2],  lw=1.0, label=pos[2] )
#axarr[1].plot(x, plot_ingredients[2], ls='-', c=set2[1],  lw=1.0, label=pos[0] )
#axarr[1].plot(x, plot_ingredients[3], ls='-', c=set2[1],  lw=1.0, label=pos[1] )


axarr[0].set_xlim(x_lim)
axarr[0].set_ylim(y_lim)
#axarr[0].set_xlabel('E - E$_F$ (eV)', fontdict=font)
axarr[0].set_ylabel(r"$E_{\mathrm{BP}}-E_{\mathrm{FM}}$/f.u. (meV)", fontdict=font)
axarr[0].legend(loc='lower left',frameon=True, prop={'size':9}, labelspacing=0.18)
axarr[0].tick_params(axis='x', colors='black',labelsize=9,width=1)
axarr[0].tick_params(axis='y', colors='black',labelsize=9,width=1)
axarr[0].grid(True)

axarr[1].set_xlim(x_lim)
axarr[1].set_ylim(y_lim2)
axarr[1].set_xlabel(r"$f_{\mathrm{SOC}}$", fontdict=font)
axarr[1].set_ylabel(r"Mag. mom. / f.u. ($\mu_{B}$)", fontdict=font)
axarr[1].legend(loc='lower left',frameon=True, prop={'size':9}, labelspacing=0.18)
axarr[1].tick_params(axis='x', colors='black',labelsize=9,width=1)
axarr[1].tick_params(axis='y', colors='black',labelsize=9,width=1)
axarr[1].grid(True)

#plt.annotate('local max', xy=(2, 1), xytext=(2.9, 1.1),
#                     arrowprops=dict(facecolor='black', shrink=0.05))

#plt.set_xlim(x_lim)
#plt.set_ylim(y_lim)
#plt.set_xlabel('E - E$_F$ (eV)', fontdict=font)
#plt.set_ylabel('LDOS (states/eV)', fontdict=font)
#plt.legend(loc='lower left',frameon=True, prop={'size':7}, labelspacing=0.18)
#plt.tick_params(axis='x', colors='black',labelsize=7,width=1)
#plt.tick_params(axis='y', colors='black',labelsize=7,width=1)



fig = matplotlib.pyplot.gcf()

fig.set_size_inches(6, 4)
filename = 'MnGe_ferro_bp.pdf' 
fig.savefig(filename,transparent=True,dpi=2540, bbox_inches='tight')
#plt.show()
