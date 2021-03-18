


# Importing modules 
print('Importing modules ...')
import numpy as np 
import matplotlib.pyplot as plt
import datetime
from datetime import timedelta
import glob
from matplotlib import gridspec
import matplotlib
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
matplotlib.rcParams.update({'font.size': 15}) 
import matplotlib.colors as colors
import pdb
import pylab


print('Loading the visibility & frequency data from numpy arrays ...')

a = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/dynamic_spectrum_from_interferometric_data_7576/visibilities.npy')
b = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/dynamic_spectrum_from_interferometric_data_7872/visibilities.npy')
a = np.delete(a, slice(586,592), axis = 3)
b = np.delete(b, slice(582,592), axis = 3)

vis  = np.concatenate((a,b), axis = 3)
  
freq = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/dynamic_spectrum_from_interferometric_data_7576/freq.npy')


print('Calculating times ...')
start_time = datetime.datetime(2014, 9, 28, 2, 46, 6, 0)
timestep = timedelta(seconds = 0.5)
times = [(start_time + i*timestep) for i in range(len(vis[0][0][0]))]
ylabels = [[131.5, 132.5, 133.5], [118.5, 119.5, 120.5], 
            [107.0, 108.0, 109.0], [97.0, 98.0, 99.0], 
            [88.0, 89.0, 90.0], [79.0, 80.0, 81.0]]
#pdb.set_trace()
fig, axs = plt.subplots(6, sharex=True, figsize=(7.4,5.8))
fig.subplots_adjust(hspace=0.02)

vmin = np.log10(1000)
vmax = np.log10(180000)
t0 = dates.date2num(times[0])
t1 = dates.date2num(times[-1])

for i in range(6):
    
    XX = vis[i][0]*np.conj(vis[i][0])
    
    YY = vis[i][3]*np.conj(vis[i][3])
    f0 = freq[i][0][0]
    f1 = freq[i][-1][0]

    spectro = np.log10(np.sqrt(np.real(XX + YY))/2)
    peak = axs[i].imshow(spectro, 
        cmap=plt.get_cmap('Spectral_r'), 
        vmin=vmin, vmax = vmax, aspect = 'auto', origin = 'lower',
        extent=(t0, t1, f0, f1))
    
    # specify number of axis label
    axs[i].set_yticks(ylabels[i])
    axs[i].set_yticklabels(ylabels[i])
    axs[i].label_outer()
    # save the plot as png


axs[5].xaxis_date()
axs[5].xaxis.set_major_locator(dates.MinuteLocator())
# specify the time format for x-axis
axs[5].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
axs[5].xaxis.set_major_locator(MaxNLocator(nbins = 8) )
axs[5].set_xlabel('Start Time: '  + str(start_time) + ' (UT)')
axs[2].set_ylabel('Frequency (MHz)')
cbar_ax = fig.add_axes([0.91, 0.1, 0.02, 0.8])
cbar = fig.colorbar(peak, cax=cbar_ax, ticks= [np.log10(1000),np.log10(2000),np.log10(5000),np.log10(10000),np.log10(20000),np.log10(40000),np.log10(80000),np.log10(160000)])
cbar.ax.set_yticklabels(['0.01','0.02','0.05','0.1','0.2','0.4','0.8','1.6'])
cbar.set_label('Jy/beam', size=10)
cbar.ax.set_title('X 10$^{5}$')
#plt.savefig(output_loc + str(start_time)[0:10]+'_dynamic_spectrum_1095907576.png') 
plt.show()
