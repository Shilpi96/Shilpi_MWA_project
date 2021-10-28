
import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import matplotlib.colors
import pdb
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
from datetime import timedelta
import glob
from matplotlib import gridspec
import matplotlib.cm as cm


### Loading the position in the x and y axis of the radio source

freq_list = ['53~55','49~50','44~47','40~43','21~23','17~20','12~15','8~11']
pos = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/new_analysis/103-104/Gauss_params.npy') 
rpos = np.zeros((len(freq_list), pos[0].shape[0], 2))
ntimes = pos[0].shape[0]
#pdb.set_trace()
for j in range(8):
	for i in range(8,pos[0].shape[0]-1):
		rpos[7-j][i][0] = pos[j][i+1][0] - pos[j][i][0] ### X shift
		rpos[7-j][i][1] = pos[j][i+1][1] - pos[j][i][1] ### Y shift

print('Calculating times ...')
start_time = datetime.datetime(2014, 9, 28, 2, 46, 6, 0)
timestep = timedelta(seconds = 1)
times = [(start_time + i*timestep) for i in range(pos[0].shape[0]-1)]
x = [dates.date2num(t) for t in times]
#X = x[skip]


###### Plotting the timeseries for shift in source position
x_data = [rpos[0][i][0] for i in range(pos[0].shape[0]-1)]
y_data = [rpos[0][i][1] for i in range(pos[0].shape[0]-1)]

x_data = np.array(x_data)
y_data = np.array(y_data)
length = np.sqrt(x_data**2 + y_data**2)
Y = [131.58,131.74,131.94,132.08,132.86,133.02,133.22,133.36]

Length = [length[i] for i in range(pos[0].shape[0]-1)]
#pdb.set_trace()
fig = plt.figure(figsize=(14, 7))
ax2 = fig.add_subplot(111)
	
ax2.plot_date(x, Length, '-', label = Y[7], color = 'Green')
ax2.xaxis_date()
ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
#ax1.set_ylabel(r'$\Delta$RA (Arcmin)')
ax2.set_ylabel(r'$\Delta$position shift (Arcsec)')
ax2.set_xlim(x[8], x[287])
ax2.set_ylim(-0.2,4)

plt.show()


