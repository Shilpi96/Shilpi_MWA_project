# Importing modules 
print('Importing modules ...')
import numpy as np 
import matplotlib.pyplot as plt
from astropy import units as u
import datetime
from datetime import timedelta
import glob
from matplotlib import gridspec
import matplotlib
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
matplotlib.rcParams.update({'font.size': 11}) 
import matplotlib.colors as colors
import pdb
import pylab
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as unit
from sunpy.coordinates.sun import sky_position as sun_position
import sunpy.coordinates.sun as sun_coord
import sunpy.map
from astropy.visualization import ImageNormalize,LogStretch
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import sunpy.visualization.colormaps as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


#### Joining the radio sources with lines on the AIA image

def joinpoints(gauss_file, fint=0, tint=0,colour = 0,color = 'Reds',linewidth = 1,a = 0):

	pos = gauss_file	

	xasec = pos[fint][tint][0]*u.arcsec
	yasec = pos[fint][tint][1]*u.arcsec
	
	xdeg = xasec.to('deg')
	ydeg = yasec.to('deg')
	
	Fint = fint
	Tint = tint+1
	
	Xasec = pos[Fint][Tint][0]*u.arcsec
	Yasec = pos[Fint][Tint][1]*u.arcsec
	
	Xdeg = Xasec.to('deg')
	Ydeg = Yasec.to('deg')
	
	xx = [xdeg,Xdeg]*u.deg
	yy = [ydeg,Ydeg]*u.deg
	
	if a == 0:
		cm = pylab.get_cmap(color)
		collist = cm(np.linspace(0, 255, 288).astype(int)) 
	
		ax1.plot(xx*u.deg, yy*u.deg,color=collist[colour],transform=ax1.get_transform("world"), linewidth = linewidth)
	else:
		ax1.plot(xx*u.deg, yy*u.deg,color=color,transform=ax1.get_transform("world"), linewidth = linewidth)
	
	return ax1




############## Defining Axes

fig = plt.figure(figsize=(9, 7))

####### Plotting the AIA map
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

ax1 = plt.subplot(111, projection=aia_171_map)
ax1.imshow(rgb)
aia_171_map.draw_grid(axes=ax1,annotate=False)
##### Setting limit in both x and y axis
xlims_world = [0, 550]*unit.arcsec
ylims_world = [-1000, -500]*unit.arcsec
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


####### Getting coordinates of the gaussian

root = '/mnt/LOFAR-PSP/shilpi_MWA_project/'

gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')

#pdb.set_trace()

for i in range(8,gauss_file1.shape[1]-80):
	joinpoints(gauss_file1, fint=0, tint=i,colour = i,color = 'black',linewidth = 2.5,a = 1)
	joinpoints(gauss_file1, fint=3, tint=i,colour = i,color = 'black',linewidth = 2.5,a = 1)

for i in range(8,gauss_file1.shape[1]-80):
	joinpoints(gauss_file1, fint=0, tint=i,colour = i,color = 'Reds',linewidth = 1.5)
	joinpoints(gauss_file1, fint=3, tint=i,colour = i,color = 'Purples', linewidth = 1.5)
	

####### Construct a colorbar
	
ax = fig.add_axes([0.72, 0.1, 0.01, 0.78])
cmap = mpl.cm.Reds
norm = mpl.colors.Normalize(vmin=1, vmax=287-70)

cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_ticks([50, 100, 150, 200])
cb1.set_ticklabels([ ])                               
#cb1.ax.tick_params(labelsize=10)
#cb1.set_label('Seconds',size = 10)

ax1 = fig.add_axes([0.75, 0.1, 0.01, 0.78])
cmap = mpl.cm.Purples
norm = mpl.colors.Normalize(vmin=1, vmax=287-70)

cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
#cb1.ax1.tick_params(labelsize=10)
cb1.set_label('Seconds',size = 10)

plt.show()


