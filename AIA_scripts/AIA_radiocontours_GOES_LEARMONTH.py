import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from datetime import datetime
from matplotlib.patches import Rectangle

import sunpy
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.time import TimeRange, parse_time
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from sunpy.coordinates import frames, sun
from astropy.visualization import ImageNormalize,LogStretch
import astropy.units as unit
from astropy.io import fits
from astropy.time import Time
import sunpy.timeseries as ts
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import pdb
import matplotlib
matplotlib.rcParams.update({'font.size': 11.5}) 


### Read the radio fits file

def read_radiofits(fits_file):
	
	hdu = fits.open(list1[i])[0]
	header = hdu.header
	data = hdu.data
	MWA_map = sunpy.map.Map(data, header)
	
	bl = SkyCoord(-3000*u.arcsec, -3000*u.arcsec, frame=MWA_map.coordinate_frame)
	tr = SkyCoord(3000*u.arcsec, 3000*u.arcsec, frame=MWA_map.coordinate_frame)
	MWA_submap = MWA_map.submap(bl, top_right=tr)
	
	return MWA_submap


####### Creating the timeseries
goes = ts.TimeSeries('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/GOES_data/sci_gxrs-l2-irrad_g15_d20140928_v0-0-0.nc')
#goes = ts.TimeSeries(files, concatenate=True) ####### if goes is returning two files for GOES 13 and 15, you want to concatenate them
goes_tr = goes.truncate("2014-09-28 02:30:00", "2014-09-28 03:40:00")

date_format_goes = dates.DateFormatter('%H:%M')

### Arranging the axes
fig = plt.figure(figsize=(9, 9))
outer = gridspec.GridSpec(3, 3, wspace=0.2, hspace=0.29)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0,:3], wspace=0.1, hspace=0.3)
inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[1,:3], wspace=0.1, hspace=0.3)
inner2 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[2,0], wspace=0.1, hspace=0.3)                
inner3 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[2,1], wspace=0.1, hspace=0.3)                
inner4 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[2,2], wspace=0.01, hspace=0.3)                



ax1 = plt.subplot(inner0[0])
gax = goes_tr.plot(ax1, legend=False, rot=0)
gax.xaxis.set_major_formatter(date_format_goes)
gax.yaxis.set_minor_locator(MultipleLocator(5))
#gax.set_xlabel("Time (UTC)",fontsize=9)
#gax.set_ylabel(r"Watts m$^{-2}$")
gax.set_xlabel("Time (UTC)")
gax.text(0.6,0.9,'GOES-15 Solar X-ray Flux', transform=gax.transAxes)
#gax.text(0.05,0.9,'A', fontdict={'size':10}, transform=gax.transAxes)
gax.set_yscale("log")
gax.tick_params(axis='x')
gax.tick_params(axis='y')
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

labels = [r"0.5-4 $\AA$",r"1-8 $\AA$"]
gax.yaxis.grid(True)
gax.legend(handles, labels, fontsize=10, loc = 'lower right')


######## Plotting Learmonth data
import os
import numpy as np
import matplotlib.dates as mdates
import datetime 
from scipy.io import readsav
from astropy.time import Time


spectro = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/learmonth_spectra.npy')
times = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/times.npy')
freqs = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/freq.npy')

#pdb.set_trace()
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


ax2 = plt.subplot(inner1[0])
ax2.pcolormesh(dtimes, freqs, spectro, 
			cmap = plt.get_cmap('Spectral_r'),
			vmin = np.percentile(spectro, 55),
			vmax = np.percentile(spectro, 99))

ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
ax2.set_xlim(t0plot, t1plot)
#ax2.set_yticklabels([])
#ax2.set_ylabel("Frequency (MHz)")
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

# Plot the frequenc and time point for the radio contours
ax2.scatter(datetime.datetime(2014, 9, 28, 2, 49, 30), freqs[621],facecolor = 'Red', edgecolor = 'black',  linewidth = 0.75, s = 60)
ax2.scatter(datetime.datetime(2014, 9, 28, 2, 49, 30), freqs[573],facecolor = 'Blue', edgecolor = 'black',  linewidth = 0.75, s = 60)
ax2.scatter(datetime.datetime(2014, 9, 28, 2, 49, 30), freqs[529],facecolor = 'Green', edgecolor = 'black',  linewidth = 0.75, s = 60)


ax2.set_aspect('auto')
ax2.text(0.8,0.06,'Learmonth', transform=ax2.transAxes, c = 'white')
#ax2.text(0.05,0.9,'B', fontdict={'size':10}, transform=ax2.transAxes, c = 'white')
#plt.savefig('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth-goes.png',dpi = 200)


######### plot the AIA image with the radio contours

root = '/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/new_analysis/'
root1 = '/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1-1095907872'
list1 = [root + '103-104/correct_fits_files/40~43/1095907576_103-104_chan_40~43_024930.0.fits', root + '093-094/correct_fits_files/40~43/1095907576_093-094_chan_40~43_024930.0.fits', root + '084-085/correct_fits_files/40~43/1095907576_084-085_chan_40~43_024930.0.fits']

MWA_submaps = []

for i in range(len(list1)):
	MWA_submaps.append(read_radiofits(list1[i]))
	

############## Defining Axes

#fig = plt.figure(figsize=(7, 9))
#fig.subplots_adjust(hspace=0.3,wspace= 0.05)


######### Getting AIA data
aia_map = sunpy.map.Map('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/aiaimages_193/aia_lev1_193a_2014_09_28t02_49_30_84z_image_lev1.fits')

for i in range(3):
	
	if i ==0: ax3 = plt.subplot(inner2[0], projection=aia_map)
	if i ==1: ax3 = plt.subplot(inner3[0], projection=aia_map)
	if i ==2: ax3 = plt.subplot(inner4[0], projection=aia_map)

	ax3.set_yticklabels([])
	cmap = ['Reds_r','Blues_r','Greens_r','Greens_r','Oranges_r','Greys_r']
	with frames.Helioprojective.assume_spherical_screen(aia_map.observer_coordinate):
		aia_map.plot(cmap = 'Greys',norm=ImageNormalize(vmin=100,vmax=4000,stretch=LogStretch(10)))
    
		
		ax3.contour(MWA_submaps[i].data, 
               levels=np.arange(0.5, 1, 0.1)*np.max(MWA_submaps[i].data),
               transform=ax3.get_transform(MWA_submaps[i].wcs), cmap=cmap[i])

##### Setting limit in both x and y axis
	xlims_world = [-100, 900]*unit.arcsec
	ylims_world = [-1200, -40]*unit.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
	pixel_coords = aia_map.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax3.set_xlim(xlims_pixel)
	ax3.set_ylim(ylims_pixel)	
	ax3.patch.set_facecolor('white')
	ax3.set_xlabel(' ')
	ax3.set_ylabel('   ')
	
	
#	if i ==0: ax3.set_ylabel(' Solar Y (arcsec)')
	if i in (0,2):ax3.set_title('  ')
	y = ax3.axes.coords[1]
	if i in (1,2):y.set_ticklabel_visible(False)
	ax3.xaxis.set_major_locator(plt.MaxNLocator(3))
#	ax3.set_yticks([-200, -400, -600, -800, -1000, -1200] *unit.arcsec)
#	ax3.set_yticklabels((" ", " "," "," "," "," "))

plt.show()

