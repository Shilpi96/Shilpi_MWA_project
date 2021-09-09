#### This code creates a plot of AIA ratio images with band-splitting sources and Learmonth (fitted curve on both bands) and MWA with overplotted sources.


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
from scipy.optimize import curve_fit
from matplotlib.patches import Rectangle
from matplotlib import cm
import matplotlib
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
	collist = cm(np.linspace(0, 255, 10).astype(int))   #### Creating the colour points in a step 10 from 0 to 255 considering you have 10 points to plot

	xdeg, ydeg = getarcsec(gauss_file, fint=fint, tint=tint)
	ax1.plot(xdeg, ydeg, marker, color=collist[colour], mec = 'Black',transform=ax1.get_transform('world'), markersize=10)
	return ax1

### Function to join the points for the same time with line

def joinpoints(gauss_HF,gauss_LF,ax1,Ufint = 0, Utint = 0,Lfint = 0, Ltint = 0,colour = 0):
	cm = pylab.get_cmap('Reds')
	collist = cm(np.linspace(0, 255, 10).astype(int))   #### Creating the colour points in a step 10 from 0 to 255 considering you have 10 points to plot

	xdeg,ydeg = getarcsec(gauss_HF, fint=Ufint, tint=Utint)
	xdeg1,ydeg1 = getarcsec(gauss_LF, fint=Lfint, tint=Ltint)
	xx = [xdeg,xdeg1]*u.deg
	yy = [ydeg,ydeg1]*u.deg
	ax1.plot(xx*u.deg, yy*u.deg,color=collist[colour],transform=ax1.get_transform("world"))
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
def plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Utpoint=0, Ufpoint=0,Ltpoint=0,Lfpoint=0,Ltpoint1=0,Lfpoint1=0,i=0,colour = 0,colour1 =0, ylabel = ' '):
	f0 = freq[i][0][0]
	f1 = freq[i][-1][0]
	
	ax = plt.Subplot(fig, inner2[i])
	peak = ax.imshow(spectra[i],cmap=plt.get_cmap('Spectral_r'),
        vmin=vmin, vmax = vmax, aspect = 'auto', origin = 'lower',
        extent=(t0, t1, f0, f1))
	
	ax.set_yticks(ylabels[i])
	ax.set_yticklabels(ylabels[i])
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
	## Plotting the points on the plot
	cm = pylab.get_cmap('Blues')
	collist = cm(np.linspace(0, 255, 10).astype(int)) 
	cm1 = pylab.get_cmap('Greens')
	collist1 = cm1(np.linspace(0, 255, 10).astype(int)) 

	if Ufpoint != 0 :ax.scatter(Utpoint,Ufpoint,facecolor=collist[colour], edgecolor='black',  linewidth = 1,s = 55)
	if Lfpoint != 0 :ax.scatter(Ltpoint,Lfpoint,facecolor=collist1[colour1], edgecolor='black', marker = 'v', linewidth = 1,s = 55)
	if Lfpoint1 != 0 :ax.scatter(Ltpoint1,Lfpoint1,facecolor=collist1[9], edgecolor='black', marker = 'v', linewidth = 1,s = 55)
	fig.add_subplot(ax)
	

### Function to overplot the points on the Learmonth
def plot_lrmthpoints(ax2,cm,cm1,Utpoint=0, Ufpoint=0,Ltpoint=0,Lfpoint=0,Ltpoint1=0,Lfpoint1=0,colour = 0,colour1 =0):
	cm = pylab.get_cmap('Blues')  ### Upper band
	collist = cm(np.linspace(0, 255, 10).astype(int)) 

	cm1 = pylab.get_cmap('Greens')  ### Lower band
	collist1 = cm1(np.linspace(0, 255, 10).astype(int)) 

	if Ufpoint != 0 :ax2.scatter(Utpoint,Ufpoint,facecolor=collist[colour], edgecolor='black',  linewidth = 1,s = 55)
	if Lfpoint != 0 :ax2.scatter(Ltpoint,Lfpoint,facecolor=collist1[colour1], edgecolor='black', marker = 'v', linewidth = 1,s = 55)
	if Lfpoint1 != 0 :ax2.scatter(Ltpoint1,Lfpoint1,facecolor=collist1[9], edgecolor='black', marker = 'v', linewidth = 1,s = 55)
	
	
######## Creating the RGB array and get the aia 171 map
print('Creating the RGB array for a specific time')
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/'

data_171 = np.load(root+'AIA_data/AIA_normailse_data/data_171.npy')
data_193 = np.load(root+'AIA_data/AIA_normailse_data/data_193.npy')
data_211 = np.load(root+'AIA_data/AIA_normailse_data/data_211.npy')

f = 83
i = 73
ratio_a = getratio(data_171, f=f,i=i)
ratio_b = getratio(data_193, f=f,i=i)
ratio_c = getratio(data_211, f=f,i=i)

rgb = np.zeros((4096, 4096, 3))
rgb[:, :, 0] = ratio_a
rgb[:, :, 1] = ratio_b
rgb[:, :, 2] = ratio_c

### load the time from the aia files
print('loading the time from the aia files')
time_171 = np.load(root+'AIA_data/AIA_normailse_data/Time_171.npy')
aia_171_map = sunpy.map.Map(root+'AIA_data/aiaimages_171/aia_lev1_171a_2014_09_28t'+time_171[f])

############## Defining Axes
print('defining the axes')
fig = plt.figure(figsize=(14, 7))
outer = gridspec.GridSpec(2, 4, wspace=0.9, hspace=0.15)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:2,:2], wspace=0.1, hspace=0.3)

##### Plot the AIA three colour ratio image
print('plotting the AIA image')
ax1 = plt.subplot(inner0[0], projection=aia_171_map)
ax1.imshow(rgb)
aia_171_map.draw_grid(axes=ax1,annotate=False)
##### Setting limit in both x and y axis
xlims_world = [-350, 1225]*unit.arcsec
ylims_world = [-1400, 200]*unit.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_171_map.coordinate_frame)
pixel_coords = aia_171_map.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value
ax1.set_xlim(xlims_pixel)
ax1.set_ylim(ylims_pixel)	
ax1.patch.set_facecolor('black')
ax1.set_xlabel('X (arcsec)')
ax1.set_ylabel('Y (arcsec)')
ax1.set_title(aia_171_map.meta['date-obs'])
##### Draw a line to do point and click
''' 
xlims_world = [354, 524]*unit.arcsec
ylims_world = [-348, -707]*unit.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_171_map.coordinate_frame)
pixel_coords = aia_171_map.world_to_pixel(world_coords)
ax1.plot(pixel_coords[0], pixel_coords[1],color="blue")'''


###########     Overplot the MWA points.
## Colour for the Upper band of the band-split
cm = pylab.get_cmap('Blues')
gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')
gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/093-094/coords_in_arcsec.npy')
gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/084-085/coords_in_arcsec.npy')
gauss_file3 = np.load(root+'CLEAN1-1095907872/new_analysis/076-077/coords_in_arcsec.npy')

overplot(gauss_file0,ax1,fint = 0, tint =185,colour = 3,color = 'Blues',marker = 'o')  ## 02:49:11
overplot(gauss_file1,ax1,fint = 0, tint =256,colour = 5,color = 'Blues',marker = 'o')  ## 02:50:23
overplot(gauss_file2,ax1,fint = 0, tint =28,colour = 7,color = 'Blues',marker = 'o')   ## 02:51:30
overplot(gauss_file3,ax1,fint = 0, tint =129,colour = 8,color = 'Blues',marker = 'o')  ## 02:53:11

## Plotting the points from the lower band of the band split
cm1 = pylab.get_cmap('Purples')
Gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/084-085/coords_in_arcsec.npy')
Gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/076-077/coords_in_arcsec.npy')
Gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/069-070/coords_in_arcsec.npy')


overplot(Gauss_file0,ax1,fint = 0, tint =185,colour = 3,color = 'Greens',marker = 'v')  ## 02:49:11
overplot(Gauss_file1,ax1,fint = 0, tint =257,colour = 5,color = 'Greens',marker = 'v')  ## 02:50:23
overplot(Gauss_file2,ax1,fint = 0, tint =32,colour = 7,color = 'Greens',marker = 'v')   ## 02:51:30
overplot(Gauss_file2,ax1,fint = 7, tint =129,colour = 9,color = 'Greens',marker = 'v')  ## 02:53:11

###### Join the MWA points with the lines
joinpoints(gauss_file0,Gauss_file0,ax1,Ufint = 0, Utint = 185,Lfint = 0, Ltint = 185,colour = 3)  ## 02:49:11
joinpoints(gauss_file1,Gauss_file1,ax1,Ufint = 0, Utint = 256,Lfint = 0, Ltint = 257,colour = 5)  ## 02:50:23
joinpoints(gauss_file2,Gauss_file2,ax1,Ufint = 0, Utint = 28,Lfint = 0, Ltint = 32,colour = 7)  ## 02:51:30
joinpoints(gauss_file3,Gauss_file2,ax1,Ufint = 0, Utint = 129,Lfint = 7, Ltint = 129,colour = 8)  ## 02:53:11

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
fpoint=[133.36, 120.56,109.04,98.80,89.84,88.06]

### Plot the MWA spectra

ax = plt.Subplot(fig, inner2[0])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Utpoint=tpoint[0], Ufpoint=fpoint[0],i = 0,colour = 3)

ax = plt.Subplot(fig, inner2[1])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Utpoint=tpoint[1], Ufpoint=fpoint[1],i = 1,colour = 5)

ax = plt.Subplot(fig, inner2[2])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Utpoint=tpoint[2], Ufpoint=fpoint[2],Ltpoint=tpoint[0],Lfpoint=fpoint[2],i = 2,colour= 7,colour1 =3, ylabel = 'Frequency (MHz) ')

ax = plt.Subplot(fig, inner2[3])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Utpoint=tpoint[3], Ufpoint=fpoint[3],Ltpoint=tpoint[1],Lfpoint=fpoint[3],i = 3,colour= 8,colour1 =5)

ax = plt.Subplot(fig, inner2[4])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,Ltpoint=tpoint[2],Lfpoint=fpoint[4],Ltpoint1=tpoint[3],Lfpoint1=fpoint[5],i = 4,colour1=7)

ax = plt.Subplot(fig, inner2[5])
plot_MWAspectra(spectra,freq,ax,ylabels,TIME,vmin,vmax,t0,t1 ,cm,cm1,i = 5)


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
##### Fit the curve on the bandsplitting on the learmonth
## Points for fitting the curve
x_coords = np.load(root +'AIA_data/npy_files/x_coord_learmonth.npy')
y_coords = np.load(root +'AIA_data/npy_files/y_coord_learmonth.npy')

## Points for fitting the curve on HFB
x_coords_HFB = np.load(root+'AIA_data/npy_files/x_coord_learmonth_HFB.npy')
y_coords_HFB = np.load(root+'AIA_data/npy_files/y_coord_learmonth_HFB.npy')

## Chosen points for HFB
x_point_HFB = np.load(root+'AIA_data/npy_files/x_coord_HFB.npy')
y_point_HFB = np.load(root+'AIA_data/npy_files/y_coord_HFB.npy')

## Chosen Points for LFB
x_point_LFB = np.load(root+'AIA_data/npy_files/x_coord_LFB.npy')
y_point_LFB = np.load(root+'AIA_data/npy_files/y_coord_LFB.npy')

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
ax2.plot(X, y_new,color = 'black' )

#### HFB
# curve fit
popt, _ = curve_fit(func, x1, y1)
# summarize the parameter values
a, b, c , d= popt
print('y = %.5f * x + %.5f * x^2 + %.5f * x^3 + %.5f'% (a, b, c, d)) 
y_new = func(x1,a,b,c,d)
#Plot the polynomial fit
ax2.plot(X1, y_new,color = 'black' )

##### Adding verticle dotted lines on the plot
for j in range(len(TIME)):
	ax.axvline(x= TIME[j], color='black', linestyle='dotted', linewidth = 0.5)

##### Overplot points on the learmonth spectra
plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[0], Ufpoint=fpoint[0],colour = 3)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[1], Ufpoint=fpoint[1],colour = 5)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[2], Ufpoint=fpoint[2],Ltpoint=tpoint[0],Lfpoint=fpoint[2],colour= 7,colour1 =3)

plot_lrmthpoints(ax2,cm,cm1,Utpoint=tpoint[3], Ufpoint=fpoint[3],Ltpoint=tpoint[1],Lfpoint=fpoint[3],colour= 8,colour1 =5,)

plot_lrmthpoints(ax2,cm,cm1,Ltpoint=tpoint[2],Lfpoint=fpoint[4],Ltpoint1=tpoint[3],Lfpoint1=fpoint[5],colour1=7)

plt.show()
