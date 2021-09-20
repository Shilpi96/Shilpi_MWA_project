#### This code creates AIA (with PFSS model) and LEARMONTH spectra with all the frequencies at same time from two bands highlighting the band splitting.


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
from scipy.optimize import curve_fit
matplotlib.rcParams.update({'font.size': 9}) 


### Function to get the positions of the radio sources
def getarcsec(gauss_file, fint=0, tint=0):

	pos = gauss_file	

	xasec = pos[fint][tint][0]*u.arcsec
	yasec = pos[fint][tint][1]*u.arcsec

#	pdb.set_trace()
	return xasec.to('deg'), yasec.to('deg')


### Function to plot the points on the AIA image
def overplot(gauss_file,ax1,fint=0, tint=0,colour = 0,color = 'Blues',marker = 'o'):

	cm = pylab.get_cmap(color)
	collist = cm(np.linspace(0, 255, 100).astype(int))   #### Creating the colour points in a step 10 from 0 to 255 considering you have 10 points to plot

	xdeg, ydeg = getarcsec(gauss_file, fint=fint, tint=tint)
	plt_points = ax1.plot(xdeg, ydeg, marker, color=collist[colour], mec = 'Black',transform=ax1.get_transform('world'), markersize=10)
	
	return ax1

### Function to plot the MWA spectra
def plot_MWAspectra(spectra, freq, ax,ylabels,TIME,vmin,vmax,t0,t1,color = 'Blues', ylabel = ' ',xlabel = ''):
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
	ax.set_xlabel(xlabel)
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
	
### Function to get the index number of the points

def get_index(radio_times, lrmnth_times, i = 8):
	delt = abs(radio_times - lrmnth_times[:, np.newaxis])
	index = []
	for i in range(i):
		delt[i] = np.where(delt[i] == delt[i].min())
		index.append(delt[i][0][0])
	return index
	
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
fig = plt.figure(figsize=(13, 6))
outer = gridspec.GridSpec(1, 2, wspace=0.15, hspace=3)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.15)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,:1], wspace=0.05, hspace=0.3)

##### Plot the AIA three colour ratio image
print('plotting the AIA image')
ax1 = plt.subplot(inner0[0], projection=aia)
aia.plot(ax1)

for fline in flines:
	color = {0: 'whitesmoke', -1: 'tab:blue', 1: 'tab:red'}.get(fline.polarity)	
	ax1.plot_coord(fline.coords, alpha=0.8, linewidth=1, color=color, zorder=1)

	
ax1.patch.set_facecolor('black')

##### Setting limit in both x and y axis
xlims_world = [-300, 1225]*unit.arcsec
ylims_world = [-1500, 50]*unit.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia.coordinate_frame)
pixel_coords = aia.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value
ax1.set_xlim(xlims_pixel)
ax1.set_ylim(ylims_pixel)
ax1.set_ylabel('Y (arcsec)')
ax1.set_xlabel('X (arcsec)')	
ax1.patch.set_facecolor('black')


###########     Overplot the MWA points.

#### Time chosen for Upper bands from the learmonth spectra
x_point_HFB = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/x_coord_HFB.npy')
HFB_time = dates.num2date(x_point_HFB)
HFB_time = [HFB_time[i].replace(tzinfo=None) for i in range(len(HFB_time))]
excluded_index = [1,19,25]
HFB_time = [i for n, i in enumerate(HFB_time) if n not in excluded_index]
HFB_time = np.array(HFB_time)
#### Time chosen for lower band from the learmonth spectra
x_point_LFB = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/x_coord_LFB.npy',allow_pickle=True)
LFB_time = dates.num2date(x_point_LFB)

LFB_time[15] = datetime.datetime(2014, 9, 28, 2, 48, 22, 655914, tzinfo=datetime.timezone.utc)
LFB_time[31] = datetime.datetime(2014, 9, 28, 2, 50, 51, 919355, tzinfo=datetime.timezone.utc)
LFB_time = [LFB_time[i].replace(tzinfo=None) for i in range(len(LFB_time))]
LFB_time = np.array(LFB_time)

## Load the time from the radio files for both the datasets
time_7576_103 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_103-104_radio_time.npy', allow_pickle=True)
time_7576_093 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_093-094_radio_time.npy', allow_pickle=True)
time_7576_084 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_084-085_radio_time.npy', allow_pickle=True)
time_7576_076 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_076-077_radio_time.npy', allow_pickle=True)
time_7872_084 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_084-085_radio_time.npy', allow_pickle=True)
time_7872_076 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_076-077_radio_time.npy', allow_pickle=True)
time_7872_069 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_069-070_radio_time.npy', allow_pickle=True)

### Getting the index of the chosen times for both bands
	
Uindex_103 = get_index(time_7576_103, HFB_time[0:7], i = 7)
Uindex_93 = get_index(time_7576_093, HFB_time[8:15], i = 7)
Uindex_84 = get_index(time_7872_084, HFB_time[18:25], i = 7)
Uindex_76 = get_index(time_7872_076, HFB_time[23:31])

Lindex_103 = get_index(time_7576_103, LFB_time[0:7], i = 7)
Lindex_93 = get_index(time_7576_093, HFB_time[8:15],i = 7)
Lindex_84 = get_index(time_7872_084, HFB_time[18:25], i = 7)



## Colour for the Upper band of the band-split
gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')
gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/093-094/coords_in_arcsec.npy')
gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/084-085/coords_in_arcsec.npy')
gauss_file3 = np.load(root+'CLEAN1-1095907872/new_analysis/076-077/coords_in_arcsec.npy')


for i in range(4):
	if i == 0:overplot(gauss_file0,ax1,fint = i, tint =Uindex_103[i],colour = i+2,color = 'Blues',marker = 'o')  ## 02:49:11
	if i >1 :overplot(gauss_file0,ax1,fint = i, tint =Uindex_103[i-1],colour = i*4,color = 'Blues',marker = 'o')  ## this is done to avoid the 49~52 fine channel of the 103-104 band.
	overplot(gauss_file1,ax1,fint = i, tint =Uindex_93[i],colour = i*4+2+20,color = 'Blues',marker = 'o')  ## 02:50:23
	overplot(gauss_file2,ax1,fint = i, tint =Uindex_84[i],colour = i*4+2+50,color = 'Blues',marker = 'o')  ## 02:51:30
#	overplot(gauss_file3,ax1,fint = i, tint =Uindex_76[i],colour = i*4+2+50,color = 'Blues',marker = 'o') ## 02:53:11
	
## Plotting the points from the lower band of the band split
Gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/084-085/coords_in_arcsec.npy')
Gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/076-077/coords_in_arcsec.npy')
Gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/069-070/coords_in_arcsec.npy')

for i in range(4):
	overplot(Gauss_file1,ax1,fint = i, tint =Lindex_93[i],colour = i*4+2+20,color = 'Greens',marker = 'v')  ## 02:50:23
	overplot(Gauss_file2,ax1,fint = i, tint =Lindex_84[i],colour = i*4+2+50,color = 'Greens',marker = 'v')  ## 02:51:30
for i in range(4):
	overplot(Gauss_file0,ax1,fint = i, tint =Uindex_103[i],colour = i*4+2,color = 'Greens',marker = 'v')	
######################
####### Defining Axes for plotting two spectras

inner2 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,1:], wspace=0.1, hspace=0.05)

##### loading the frequency and time
##### Loading the x_coords and y_coords
## Points for fitting the curve on LFB
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data'

x_coords = np.load(root +'/npy_files/x_coord_learmonth.npy')
y_coords = np.load(root +'/npy_files/y_coord_learmonth.npy')

## Points for fitting the curve on HFB
x_coords_HFB = np.load(root+'/npy_files/x_coord_learmonth_HFB.npy')
y_coords_HFB = np.load(root+'/npy_files/y_coord_learmonth_HFB.npy')

## Chosen points for HFB
x_point_HFB = np.load(root+'/npy_files/x_coord_HFB.npy')
y_point_HFB = np.load(root+'/npy_files/y_coord_HFB.npy')

## Chosen Points for LFB
x_point_LFB = np.load(root+'/npy_files/x_coord_LFB.npy')
y_point_LFB = np.load(root+'/npy_files/y_coord_LFB.npy')

## Loading the spectra data
spectro = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/learmonth_spectra.npy')
times = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/times.npy')
freqs = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/GOES_learmonth/learmonth_data/freq.npy')


#The Following takes care of the difference in epochs for time in IDL and Python.
idl_ints_epoch = datetime.datetime(1979, 1, 1, 0, 0, 0)  
unix_tai_epoch = datetime.datetime(1970, 1, 1, 1, 0, 0)
tdelt = (idl_ints_epoch - unix_tai_epoch).total_seconds()
times = [t+tdelt for t in times]
times = np.array(times)
dtimes = [datetime.datetime.fromtimestamp(t) for t in times]
start = dates.date2num(dtimes[5492])
end = dates.date2num(dtimes[5686])
width = end - start


# Choose plot times
t0plot = datetime.datetime(2014, 9, 28, 2, 46, 6)
t1plot = datetime.datetime(2014, 9, 28, 2, 55, 48)


ax = plt.subplot(inner2[0])
ax.pcolormesh(dtimes, freqs, spectro, 
			cmap = plt.get_cmap('Spectral_r'),
			vmin = np.percentile(spectro, 55),
			vmax = np.percentile(spectro, 99))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
ax.set_xlim(t0plot, t1plot)
ax.set_ylim(freqs[415], freqs[629])
#ax.scatter(x_point_HFB,y_point_HFB,color = 'blue', edgecolor='black',linewidth = 0.5,s = [30 for i in range(len(x_point_HFB))])
#ax.scatter(x_point_LFB,y_point_LFB,color = 'green', edgecolor='black',marker = 'v',s = [30 for i in range(len(x_point_LFB))])
ax.xaxis_date()
ax.xaxis.set_major_locator(dates.MinuteLocator())
# specify the time format for x-axis
ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
ax.set_xlabel('Time (UT)')
ax.set_ylabel('Frequency (MHz)')
##### Drawing transparent rectangles to show the MWA channels on the Learmonth spectra
freq = [131.2,118.4,106.88,97.0,87.68,78.72]

for i in range(len(freq)):
	rect = Rectangle((start, freq[i]), width, 2.52, facecolor='white', alpha = 0.3, edgecolor='black', linewidth = 1)
	ax.add_patch(rect) 

## Plot the points for both bands
## Plot the points for both bands
cm = pylab.get_cmap('Blues')  ### Upper band
collist = cm(np.linspace(0, 255, 100).astype(int)) 
cm1 = pylab.get_cmap('Greens')  ### Lower band
collist1 = cm1(np.linspace(0, 255, 100).astype(int)) 
for j in range(4):
	ax.scatter(x_point_HFB[0:7][j],y_point_HFB[0:7][j],facecolor=collist[j*4+2], edgecolor='black',  linewidth = 0.5, s = 65)
#	ax.scatter(x_point_LFB[0:7][j],y_point_LFB[0:7][j],facecolor=collist1[j*4+5], edgecolor='black', marker = 'v', linewidth = 0.5, s = 50)
for j in range(4):
	ax.scatter(x_point_HFB[8:15][j],y_point_HFB[8:15][j],facecolor=collist[j*4+2+20], edgecolor='black',  linewidth = 0.5, s = 65)
#	ax.scatter(x_point_LFB[8:16][j],y_point_LFB[8:16][j],facecolor=collist1[j*4+2], edgecolor='black', marker = 'v', linewidth = 0.5, s = 50)
for j in range(4):
	ax.scatter(x_point_HFB[0:7][j],y_point_LFB[18:25][j],facecolor=collist1[j*4+2], edgecolor='black', marker = 'v', linewidth = 1, s = 65)
	ax.scatter(x_point_HFB[15:24][j],y_point_HFB[15:24][j],facecolor=collist[j*4+2+50], edgecolor='black', linewidth = 1, s = 65)

for j in range(4):
	ax.scatter(x_point_HFB[8:15][j],y_point_LFB[24:31][j],facecolor=collist1[j*4+2+20], edgecolor='black', marker = 'v', linewidth = 1, s = 65)
		
#	ax.scatter(x_point_HFB[24:32][j],y_point_HFB[24:32][j],facecolor=collist[j*4+2+50], edgecolor='black', linewidth = 1, s = 50)

for j in range(4):
	ax.scatter(x_point_HFB[18:25][j],y_point_LFB[32:40][j],facecolor=collist1[j*4+2+50], edgecolor='black', marker = 'v', linewidth = 1, s = 65)

# Draw a curve fitting the points
## LFB
X = x_coords
x = X -X[0] 
y = y_coords

## HFB
X1 = x_coords_HFB
x1 = X1-X1[0] 
y1 = y_coords_HFB

#pdb.set_trace()
##### define the true objective function
def func(x, a, b, c, d):
    return a * x + b * x**2 + c * x**3 + d #* x**4 + e
#### LFB
# curve fit
popt, _ = curve_fit(func, x, y)
# summarize the parameter values
a, b, c , d= popt
print('y = %.5f * x + %.5f * x^2 + %.5f * x^3 + %.5f'% (a, b, c, d))  # * x^4 + %.5f' % (a, b, c, d, e))

y_new = func(x,a,b,c,d)

#Plot the polynomial fit
ax.plot(X, y_new,color = 'black' )


#### HFB
# curve fit
popt, _ = curve_fit(func, x1, y1)
# summarize the parameter values
a, b, c , d= popt
print('y = %.5f * x + %.5f * x^2 + %.5f * x^3 + %.5f'% (a, b, c, d)) 

y_new = func(x1,a,b,c,d)

#Plot the polynomial fit
ax.plot(X1, y_new,color = 'black' )


plt.show()
