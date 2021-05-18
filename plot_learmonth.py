

import os
import numpy as np
import matplotlib.pyplot as plt
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


file = '/Users/shilpibhunia/Documents/MWA_Project/Data/LEA_20140927_2211_srs.sav'

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
t1plot = datetime.datetime(2014, 9, 28, 3, 0)

spectro = backsub(spectro)

fig, ax = plt.subplots(1, 1)
ax.pcolormesh(dtimes, freqs, spectro, 
			cmap = plt.get_cmap('Spectral_r'),
			vmin = np.percentile(spectro, 55),
			vmax = np.percentile(spectro, 99))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
ax.set_xlim(t0plot, t1plot)
plt.show()
