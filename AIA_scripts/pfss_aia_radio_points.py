#### This code creates AIA and both spectras with all the frequencies at same time from two bands highlighting the band splitting


import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib.gridspec as gridspec
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
import pdb
from matplotlib import pylab
import matplotlib.colors as colors
import sunpy.map
from sunpy.coordinates import frames, sun
from scipy import optimize
from astropy.modeling.functional_models import Gaussian2D
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as unit
from astropy.visualization import ImageNormalize,LogStretch
import pylab
from numpy import *    ### it is needed to define ravel
import datetime
from datetime import timedelta
from matplotlib import dates
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import os
from scipy.io import readsav
from matplotlib.patches import Rectangle
from matplotlib import cm
import matplotlib
import pfsspy
import pfsspy.tracing as tracing
import astropy.constants as const
matplotlib.rcParams.update({'font.size': 9}) 


### Function to get the positions of the radio sources
def getarcsec(gauss_file, fint=0, tint=0):

	pos = gauss_file	

	xasec = pos[fint][tint][0]*u.arcsec
	yasec = pos[fint][tint][1]*u.arcsec

	return xasec.to('deg'), yasec.to('deg')


### Function to plot the points on the AIA image
def overplot(gauss_file,ax1,fint=0, tint=0,colour = 0,color = 'Blues',marker = 'o'):

	cm = pylab.get_cmap(color)
	collist = cm(np.linspace(0, 255, 100).astype(int))   #### Creating the colour points in a step 10 from 0 to 255 considering you have 10 points to plot

	xdeg, ydeg = getarcsec(gauss_file, fint=fint, tint=tint)
	plt_points = ax1.plot(xdeg, ydeg, marker, color=collist[colour], mec = 'Black',transform=ax1.get_transform('world'))
	
	return ax1

### function to get the ratio data
def getratio(data, f=0,i=0, vmin=0.5, vmax=1.5):

	
        #----------------------#
        #     Get ratio
        #
	ratio = data[f]/data[i]

        #------------------------------------#
        #     Remove NaNs and infinites
        #
	ratio = np.nan_to_num(ratio)
	ratio[np.where(ratio>1e3)] = 1.0
	ratio[np.where(ratio<-1e3)] = 1.0
	ratio = np.clip(ratio, vmin, vmax)
        
        #------------------------------------#
        #     Map the intensity scale to between 0 -> 1 for the RGB array
        #     (The RGB plotter expects intensity values from 0 to 1.
	imscale = np.linspace(vmin, vmax, len(ratio))
	newscale = np.linspace(0.0, 1.0, len(ratio))
	iscale_img = np.interp(ratio, imscale, newscale)
	return iscale_img

### Function to plot the MWA spectra
def plot_MWAspectra(spectra, freq, ax,ylabels,TIME,vmin,vmax,t0,t1,color = 'Blues', ylabel = ' '):
	f0 = freq[0][0]
	f1 = freq[-1][0]
	
	#ax = plt.Subplot(fig, inner2[i])
	peak = ax.imshow(spectra,cmap=plt.get_cmap('Spectral_r'),
        vmin=vmin, vmax = vmax, aspect = 'auto', origin = 'lower',
        extent=(t0, t1, f0, f1))
	
	ax.set_yticks(ylabels)
	ax.set_yticklabels(ylabels)
	ax.set_xticklabels([])
	## Adding verticle lines on the spectra
	for j in range(len(TIME)):
		ax.axvline(x= TIME[j], color='black', linestyle='dotted', linewidth = 0.5)
	
		
	ax.xaxis_date()
	ax.xaxis.set_major_locator(dates.MinuteLocator())
	## specify the time format for x-axis
	ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
	ax.set_xlabel('Time (UT)')
	ax.set_ylabel(ylabel)
	
	## Making the colormaps to overplot the MWA points on the spectra
	cm = pylab.get_cmap(color)  ### Upper band
	collist = cm(np.linspace(0, 255, 20).astype(int)) 

	fig.add_subplot(ax)
	

### Function to overplot the points on the Learmonth
def plot_lrmthpoints(ax2,cm,cm1,Utpoint=0, Ufpoint=0,Ltpoint=0,Lfpoint=0,Ltpoint1=0,Lfpoint1=0,colour = 0,colour1 =0):
	cm = pylab.get_cmap('Blues')  ### Upper band
	collist = cm(np.linspace(0, 255, 10).astype(int)) 

	cm1 = pylab.get_cmap('Greens')  ### Lower band
	collist1 = cm1(np.linspace(0, 255, 10).astype(int)) 

	if Ufpoint != 0 :ax2.scatter(Utpoint,Ufpoint,facecolor=collist[colour], edgecolor='black',  linewidth = 1)
	if Lfpoint != 0 :ax2.scatter(Ltpoint,Lfpoint,facecolor=collist1[colour1], edgecolor='black', marker = 'v', linewidth = 1)
	if Lfpoint1 != 0 :ax2.scatter(Ltpoint1,Lfpoint1,facecolor=collist1[9], edgecolor='black', marker = 'v', linewidth = 1)
	
	
######## Plotting the magnetic field on AIA plot

### Load a GONG magnetic field map
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/'
gong_fname = root+'pfss_gong/mrbqs140928t0214c2155_163.fits'

### We can now use SunPy to load the GONG fits file, and extract the magnetic field data.The mean is subtracted to enforce div(B) = 0 on the solar surface: n.b. it is not obvious this is the correct way to do this, so use the following lines at your own risk!

gong_map = sunpy.map.Map(gong_fname)
# Remove the mean
gong_map = sunpy.map.Map(gong_map.data - np.mean(gong_map.data), gong_map.meta)

### The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where rho = ln(r), and r is the standard spherical radial coordinate. We need to define the number of grid points in rho, and the source surface radius.

nrho = 25
rss = 2.5

### From the boundary condition, number of radial grid points, and source surface, we now construct an Input object that stores this information
input = pfsspy.Input(gong_map, nrho, rss)

## Now we construct a 5 x 5 grid of footpoitns to trace some magnetic field lines from. These coordinates are defined in the native Carrington coordinates of the input magnetogram.
# Create 5 points spaced between sin(lat)={0.35, 0.55}
s = np.linspace(-0.6, -0.1, 13)   ## (-1, 1)
# Create 5 points spaced between long={211, 262} degrees
phi = np.linspace(218, 263, 13)
print(f's = {s}')
print(f'phi = {phi}')
# Make a 2D grid from these 1D points
s, phi = np.meshgrid(s, phi)

# Now convert the points to a coordinate object
lat = np.arcsin(s) * u.rad
lon = phi * u.deg
seeds = SkyCoord(lon.ravel(), lat.ravel(), 1.01 * const.R_sun,
                 frame=gong_map.coordinate_frame)

### Compute the PFSS solution from the GONG magnetic field input
output = pfsspy.pfss(input)

### Trace field lines from the footpoints defined above.
tracer = tracing.PythonTracer()
flines = tracer.trace(seeds, output)

####### Plotting the magnetic field lines on the AIA map
aia = sunpy.map.Map('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/aiaimages_193/aia_lev1_193a_2014_09_28t02_49_30_84z_image_lev1.fits')

############## Defining Axes
print('defining the axes')
fig = plt.figure(figsize=(14, 7))
outer = gridspec.GridSpec(2, 4, wspace=0.9, hspace=0.15)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:2,:2], wspace=0.1, hspace=0.3)

##### Plot the AIA three colour ratio image
print('plotting the AIA image')
ax1 = plt.subplot(inner0[0], projection=aia)
aia.plot(ax1)

for fline in flines:
	color = {0: 'whitesmoke', -1: 'tab:blue', 1: 'tab:red'}.get(fline.polarity)	
	ax1.plot_coord(fline.coords, alpha=0.8, linewidth=1, color=color, zorder=1)

	
ax1.patch.set_facecolor('black')

##### Setting limit in both x and y axis
xlims_world = [-350, 1225]*unit.arcsec
ylims_world = [-1700, 200]*unit.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia.coordinate_frame)
pixel_coords = aia.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value
ax1.set_xlim(xlims_pixel)
ax1.set_ylim(ylims_pixel)	
ax1.patch.set_facecolor('black')

#ax1.set_title(aia.meta['date-obs'])


###########     Overplot the MWA points.
## Colour for the Upper band of the band-split
gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')
gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/093-094/coords_in_arcsec.npy')
gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/084-085/coords_in_arcsec.npy')
gauss_file3 = np.load(root+'CLEAN1-1095907872/new_analysis/076-077/coords_in_arcsec.npy')

#Ugauss_file = [gauss_file0, gauss_file1, gauss_file2, gauss_file3]


for i in range(8):
	if i == 0:overplot(gauss_file0,ax1,fint = i, tint =185,colour = i+2,color = 'Blues',marker = 'o')  ## 02:49:11
	if i >1 :overplot(gauss_file0,ax1,fint = i, tint =185,colour = i*2,color = 'Blues',marker = 'o')  ## this is done to avoid the 49~52 fine channel of the 103-104 band.
	overplot(gauss_file1,ax1,fint = i, tint =256,colour = i*2+2+14,color = 'Blues',marker = 'o')  ## 02:50:23
	overplot(gauss_file2,ax1,fint = i, tint =28,colour = i*2+2+30,color = 'Blues',marker = 'o')  ## 02:51:30
	overplot(gauss_file3,ax1,fint = i, tint =129,colour = i*2+2+50,color = 'Blues',marker = 'o') ## 02:53:11
## Plotting the points from the lower band of the band split
Gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/084-085/coords_in_arcsec.npy')
Gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/076-077/coords_in_arcsec.npy')
Gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/069-070/coords_in_arcsec.npy')

for i in range(8):
	overplot(Gauss_file0,ax1,fint = i, tint =185,colour = i*2+2,color = 'Greens',marker = 'v')  ## 02:49:11
	overplot(Gauss_file1,ax1,fint = i, tint =257,colour = i*2+2+14,color = 'Greens',marker = 'v')  ## 02:50:23
	overplot(Gauss_file2,ax1,fint = i, tint =32,colour = i*2+2+30,color = 'Greens',marker = 'v')  ## 02:51:30
	overplot(Gauss_file2,ax1,fint = i, tint =129,colour = i*2+2+50,color = 'Greens',marker = 'v')  ## 02:53:11
	
######################
####### Defining Axes for plotting two spectras
inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,2:], wspace=0.1, hspace=0.04)           ##### from 0th row and 1st column  

inner2 = gridspec.GridSpecFromSubplotSpec(6, 1,
        subplot_spec=outer[1:2,2:], wspace=0.1, hspace=0.05)


######## plotting the MWA dynamic spectra
spectra = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/1095907576/dynamic_spectrum_from_interferometric_data/spectra.npy')
  
freq = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/1095907576/dynamic_spectrum_from_interferometric_data/freq.npy')

print('Calculating times ...')
start_time = datetime.datetime(2014, 9, 28, 2, 46, 6, 0)
timestep = timedelta(seconds = 0.5)
length  = spectra.shape[2]
times = [(start_time + i*timestep) for i in range(length)]
times_mpl = [dates.date2num(t) for t in times]
       
ylabels = [[132.5], [119.5], 
            [108.0], [ 98.0], 
            [ 89.0], [ 80.0]]

vmin = np.log10(1000)
vmax = np.log10(180000)
t0 = dates.date2num(times[0])
t1 = dates.date2num(times[-1])
TIME = [(datetime.datetime(2014, 9, 28, 2, 47) +i*timedelta(seconds = 60)) for i in range(9)]
tpoint = [datetime.datetime(2014, 9, 28, 2, 49,11),datetime.datetime(2014, 9, 28, 2, 50,23),datetime.datetime(2014, 9, 28, 2, 51,30),datetime.datetime(2014, 9, 28, 2, 53,11)]
Ufpoint=[133.36, 133.02, 132.86, 132.08, 131.94, 131.74, 131.58]
Ufpoint1=[120.56,120.42,120,22,120.06,119.28,119.14,118.94,118.78]
Lfpoint=[109.04,108.9,108.7,108.54,107.76,107.62,107,42,107.26]
Lfpoint1=[98.8,98.66,98.46,98.3,97.52,97.38,97.18,97.02]
Lfpoint2=[89.84,89.7,89.5,89.34,89.56,88.42,88.22,88.06]
### Plot the MWA spectra
cm = pylab.get_cmap('Blues')  ### Upper band
collist = cm(np.linspace(0, 255, 100).astype(int)) 
cm1 = pylab.get_cmap('Greens')  ### Lower band
collist1 = cm1(np.linspace(0, 255, 100).astype(int)) 


for i in range(6):
	ax = plt.Subplot(fig, inner2[i])
	if i == 0: 
		for j in range(7):
			ax.scatter(tpoint[0],Ufpoint[j],facecolor=collist[j*2+2], edgecolor='black',  linewidth = 1)
	if i == 1: 
		for j in range(8):
			ax.scatter(tpoint[1],Ufpoint1[j],facecolor=collist[j*2+2+14], edgecolor='black',  linewidth = 1)
			
	if i == 2: 
		for j in range(8):
			ax.scatter(tpoint[0],Lfpoint[j],facecolor=collist1[j*2+2], edgecolor='black', marker = 'v', linewidth = 1)
			ax.scatter(tpoint[2],Lfpoint[j],facecolor=collist[j*2+2+30], edgecolor='black', linewidth = 1)
	if i == 3: 
		for j in range(8):
			ax.scatter(tpoint[1],Lfpoint1[j],facecolor=collist1[j*2+2+14], edgecolor='black', marker = 'v', linewidth = 1)
		for j in range(8):
			ax.scatter(tpoint[3],Lfpoint1[j],facecolor=collist[j*2+2+50], edgecolor='black', linewidth = 1)
			
	if i == 4: 
		for j in range(8):
			ax.scatter(tpoint[2],Lfpoint2[j],facecolor=collist1[j*2+2+30], edgecolor='black', marker = 'v', linewidth = 1)
		for j in range(8):
			ax.scatter(tpoint[3],Lfpoint2[j],facecolor=collist1[j*2+2+50], edgecolor='black', marker = 'v', linewidth = 1)
	
	plot_MWAspectra(spectra[i],freq[i],ax,ylabels[i],TIME,vmin,vmax,t0,t1)


######## Plotting the Learmonth dynamic spectra

spectro = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/learmonth_spectra.npy')
times1 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/times.npy')
freqs = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/freq.npy')

#The Following takes care of the difference in epochs for time in IDL and Python.
idl_ints_epoch = datetime.datetime(1979, 1, 1, 0, 0, 0)  
unix_tai_epoch = datetime.datetime(1970, 1, 1, 1, 0, 0)
tdelt = (idl_ints_epoch - unix_tai_epoch).total_seconds()
times1 = [t+tdelt for t in times1]
times1 = np.array(times1)

dtimes = [datetime.datetime.fromtimestamp(t) for t in times1]
start = dates.date2num(dtimes[5492])
end = dates.date2num(dtimes[5686])
width = end - start
# Choose plot times
t0plot = datetime.datetime(2014, 9, 28, 2, 46, 6)
t1plot = datetime.datetime(2014, 9, 28, 2, 55,48)


ax2 = plt.subplot(inner1[0])
ax2.pcolormesh(dtimes, freqs, spectro, cmap = plt.get_cmap('Spectral_r'),vmin = np.percentile(spectro, 80),vmax = np.percentile(spectro, 99))

ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#ax2.xaxis.tick_top()
ax2.set_xlim(t0plot, t1plot)
ax2.set_ylim(freqs[415], freqs[629])
ax2.set_xlabel('Time (UT)')
ax2.set_ylabel('Frequency (MHz)')

##### Drawing transparent rectangles to show the MWA channels on the Learmonth spectra
freq = [131.2,118.4,106.88,97.0,87.68,78.72]
for i in range(len(freq)):
	rect = Rectangle((start, freq[i]), width, 2.52, facecolor='white', alpha = 0.3, edgecolor='black', linewidth = 1)
	ax2.add_patch(rect) 

##### Adding verticle dotted lines on the plot
for j in range(len(TIME)):
	ax.axvline(x= TIME[j], color='black', linestyle='dotted', linewidth = 0.5)

##### Overplot points on the learmonth spectra
'''
plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[0], Ufpoint=fpoint[0],colour = 3)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[1], Ufpoint=fpoint[1],colour = 5)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[2], Ufpoint=fpoint[2],Ltpoint=tpoint[0],Lfpoint=fpoint[2],colour= 7,colour1 =3)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[3], Ufpoint=fpoint[3],Ltpoint=tpoint[1],Lfpoint=fpoint[3],colour= 8,colour1 =5,)

plot_lrmthpoints(ax2,cm,cm1,Ltpoint=tpoint[2],Lfpoint=fpoint[4],Ltpoint1=tpoint[3],Lfpoint1=fpoint[5],colour1=7)'''

plt.show()
