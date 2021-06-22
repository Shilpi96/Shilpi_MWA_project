import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from datetime import datetime
from matplotlib.patches import Rectangle

import sunpy
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.time import TimeRange, parse_time
import sunpy.timeseries as ts
import pdb

####### Downloading the GOES Data

'''tr = TimeRange(['2014-09-28 02:00', '2014-09-28 04:00'])
results = Fido.search(a.Time(tr), a.Instrument.xrs, a.goes.SatelliteNumber(15))
files = Fido.fetch(results, path='/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/GOES_data', overwrite=False)
print(results)'''

####### Creating the timeseries
goes = ts.TimeSeries('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/GOES_data/sci_gxrs-l2-irrad_g15_d20140928_v0-0-0.nc')
#goes = ts.TimeSeries(files, concatenate=True) ####### if goes is returning two files for GOES 13 and 15, you want to concatenate them
goes_tr = goes.truncate("2014-09-28 02:30:00", "2014-09-28 03:40:00")

date_format_goes = dates.DateFormatter('%H:%M')

fig = plt.figure(figsize=(6,10))
fig.subplots_adjust(hspace=0.2)
ax1 = fig.add_subplot(2,1,1)
gax = goes_tr.plot(ax1, legend=False, fontsize=16, rot=0)
 

gax.xaxis.set_major_formatter(date_format_goes)
gax.yaxis.set_minor_locator(MultipleLocator(5))
#gax.set_xlabel("Time (UTC)",fontsize=9)
gax.set_ylabel(r"Watts m$^{-2}$",fontsize=10)
gax.set_xlabel("Time (UTC)",fontsize=10)
gax.text(0.6,0.9,'GOES-15 Solar X-ray Flux', fontdict={'size':10}, transform=gax.transAxes)
gax.text(0.05,0.9,'A', fontdict={'size':10}, transform=gax.transAxes)
gax.set_yscale("log")
gax.tick_params(axis='x', labelsize=10)
gax.tick_params(axis='y', labelsize=10)
gax.autoscale(axis="x", tight=True)
gax.set_aspect('auto')
gax.axvline(x = goes_tr.index[293], color = 'black',linestyle='--')
gax.axvline(x = goes_tr.index[879], color = 'black',linestyle='--')


gax2 = ax1.twinx()
gax2.set_yscale("log")
gax.set_ylim(1e-9, 1e-3)
gax2.set_ylim(1e-9, 1e-3)
gax2.set_yticks((1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3))
gax2.set_yticklabels(('A', 'B', 'C', 'M', 'X', ' '),fontsize=10 )
gax2.set_aspect('auto')

handles, labels = gax.get_legend_handles_labels()
#handles.reverse()

labels = [r"1-8 $\AA$", r"0.5-4 $\AA$"]
gax.yaxis.grid(True)
gax.legend(handles, labels, fontsize=10, loc = 'lower right')


######## Plotting Learmonth data
import os
import numpy as np
import matplotlib.dates as mdates
import datetime 
from scipy.io import readsav
from astropy.time import Time


def backsub(data, percentile=5.0):

    # Get time slices with standard devs in the bottom nth percentile.
    # Get average spectra from these time slices.
    # Devide through by this average spec.
    # Expects (row, column)

    print('Performing background subtraction.')
    data = np.log10(data)
    data[np.where(np.isinf(data)==True)] = 0.0
    data_std = np.std(data, axis=0)
    data_std = data_std[np.nonzero(data_std)]
    min_std_indices = np.where( data_std < np.percentile(data_std, percentile) )[0]
    min_std_spec = data[:, min_std_indices]
    min_std_spec = np.mean(min_std_spec, axis=1)
    data = np.transpose(np.divide( np.transpose(data), min_std_spec))
    print('Background subtraction finished.')

    #Alternative: Normalizing frequency channel responses using median of values.
        #for sb in np.arange(data.shape[0]):
        #       data[sb, :] = data[sb, :]/np.mean(data[sb, :])

    return data


file = '/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/LEA_20140927_2211_srs.sav'

sav_data = readsav(file)
spectro = sav_data['rstn_struct']['data'][0]
times = sav_data['rstn_struct']['times'][0]
freqs = sav_data['rstn_struct']['freqs'][0]

#The Following takes care of the difference in epochs for time in IDL and Python.
idl_ints_epoch = datetime.datetime(1979, 1, 1, 0, 0, 0)  
unix_tai_epoch = datetime.datetime(1970, 1, 1, 1, 0, 0)
tdelt = (idl_ints_epoch - unix_tai_epoch).total_seconds()
times = [t+tdelt for t in times]
times = np.array(times)
dtimes = [datetime.datetime.fromtimestamp(t) for t in times]

# Choose plot times
t0plot = datetime.datetime(2014, 9, 28, 2, 40)
t1plot = datetime.datetime(2014, 9, 28, 3, 00)

spectro = backsub(spectro)

ax2 = fig.add_subplot(2,1,2)
ax2.pcolormesh(dtimes, freqs, spectro, 
			cmap = plt.get_cmap('Spectral_r'),
			vmin = np.percentile(spectro, 55),
			vmax = np.percentile(spectro, 99))

ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
ax2.set_xlim(t0plot, t1plot)
ax2.set_xlabel("Time (UTC)",fontsize=10)
ax2.set_ylabel("Frequency (MHz)",fontsize=10)
# convert to matplotlib date representation
start = dates.date2num(dtimes[5492])
end = dates.date2num(dtimes[5686])
width = end - start

# Plot rectangle
rect = Rectangle((start, 131.5), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
rect1 = Rectangle((start, 118.5), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
rect2 = Rectangle((start, 107.0), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
rect3 = Rectangle((start, 97.0), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
rect4 = Rectangle((start, 88.0), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
rect5 = Rectangle((start, 79.0), width, 2, linewidth = 0.6,edgecolor='black', facecolor='none')
ax2.add_patch(rect) 
ax2.add_patch(rect1)
ax2.add_patch(rect2)
ax2.add_patch(rect3)
ax2.add_patch(rect4)
ax2.add_patch(rect5)
ax2.set_aspect('auto')
ax2.text(0.8,0.06,'Learmonth', fontdict={'size':10}, transform=ax2.transAxes)
ax2.text(0.05,0.9,'B', fontdict={'size':10}, transform=ax2.transAxes)
plt.savefig('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth-goes.png',dpi = 200)
plt.show()

