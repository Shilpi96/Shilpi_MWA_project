##### This code creates consecutive AIA ratio images with 171,193 and 211 channels.

import glob
import pdb
import sunpy.map
from aiapy.calibrate import register, update_pointing, normalize_exposure
import os
from astropy.convolution import convolve
from astropy.convolution import Box2DKernel
import astropy.units as unit
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
from datetime import datetime
from astropy.io import fits
import matplotlib.pyplot as plt

######## building the functions
### Get the exposure time of the file
def exptime(fits_file):
	hdu = fits.open(fits_file)[1]
	header = hdu.header
	exptime = header['exptime']
	return exptime
	
### Get the dateobs of the file
def gettime(fits_file):
	hdu = fits.open(fits_file)[1]
	header = hdu.header
	TIME = header['DATE-OBS']
	time =  Time(TIME)
	time = datetime.strptime(time.value[:-1], "%Y-%m-%dT%H:%M:%S.%f")
	time = datetime.timestamp(time)
	return time


def normalise(aiamap):
	#--------------------------------#
	#    Normalise the arrays by dividing the data by the exp time to get DN/s. Very important because each AIA fits file has different exposure time.
	#    Then we clip the data by mentioning a range so that we are not working with any saturated data.
	
	data = np.array(aiamap.data/aiamap.exposure_time)

	return data

def getratio(img, imgpre, vmin=0.5, vmax=1.5):

	
        #----------------------#
        #     Get ratio
        #
	ratio = img/imgpre

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

### Get the box car smooth data
def smooth(fits_file):
	hdu = fits.open(fits_file)[1]
	kernel = Box2DKernel(10)
	data = convolve(hdu.data, kernel)
	return data

######## Getting the files
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/'
aia_211 = np.array(sorted(glob.glob(root+'aiaimages_211/*.fits')))
aia_171 = np.array(sorted(glob.glob(root+'aiaimages_171/*.fits')))
aia_193 = np.array(sorted(glob.glob(root+'aiaimages_193/*.fits')))

######## Check the images to make sure we're not using AEC-affected images and delete these AEC-affected images from the list
min_exp_t_193 = 1.0
min_exp_t_211 = 1.5
min_exp_t_171 = 1.5

## 211 channel
for i in range(len(aia_211)):
	Exptime = exptime(aia_211[i]) 
	if aia.meta['exptime']<Exptime : 
		aia_211.remove(aia_211[i])
	else: 
		print('there is no AEC effected images in 211')

## 193 channel
for i in range(len(aia_193)):
	Exptime = exptime(aia_193[i])
	if aia.meta['exptime']<Exptime : 
		aia_193.remove(aia_193[i])
	else: 
		print('there is no AEC effected images in 193')

## 171 channel
for i in range(len(aia_171)):
	Exptime = exptime(aia_171[i])
	if aia.meta['exptime']<Exptime : 
		aia_171.remove(aia_171[i])
	else: 
		print('there is no AEC effected images in 171')


######## Get the images with similar times in each channel


time211 = [gettime(f) for f in aia_211]
time211 = np.array(time211)
time193 = [gettime(f) for f in aia_193]
time193 = np.array(time193)
time171 = [gettime(f) for f in aia_171]
time171 = np.array(time171)

######## Fine the channel which has least number of files
a = [time211,time171,time193]
b = [len(a[i]) for i in range(len(a))]
print(b.index(min(b)))
#print(a.index(min(a)))

######## Create the list of files with same time for each of the channel
print('Creating the list of files with the same time')
indices_211 = np.zeros(len(time193))
indices_171 = np.zeros(len(time193))

for index, tim in enumerate(time193):
	delt211 = abs(time211 - tim)
	delt171 = abs(time171 - tim)
	
	closest211index = np.where(delt211 == delt211.min())
	closest171index = np.where(delt171 == delt171.min())

	indices_211[index] = closest211index[0][0]
	indices_171[index] = closest171index[0][0]
	
indices_211 = indices_211.astype('int')
indices_171 = indices_171.astype('int')

aia_211 = aia_211[indices_211]
aia_171 = aia_171[indices_171]

##### Creating RGB array of ration data for each channel				
####### Creating a loop to get the ratio image for all times
print('Creating the RGB array')

for i in range(len(time193)):

	aia_171_map_pre = sunpy.map.Map(aia_171[i])
	aia_193_map_pre = sunpy.map.Map(aia_193[i])
	aia_211_map_pre = sunpy.map.Map(aia_211[i])

	aia_171_map = sunpy.map.Map(aia_171[i+10])    #### taking every 10th file to get the ratio
	aia_193_map = sunpy.map.Map(aia_193[i+10])
	aia_211_map = sunpy.map.Map(aia_211[i+10])

	data_a_pre = normalise(aia_171_map_pre)
	data_b_pre = normalise(aia_193_map_pre)
	data_c_pre = normalise(aia_211_map_pre)

	data_a = normalise(aia_171_map)
	data_b = normalise(aia_193_map)
	data_c = normalise(aia_211_map)

	ratio_a = getratio(data_a, data_a_pre)
	ratio_b = getratio(data_b, data_b_pre)
	ratio_c = getratio(data_c, data_c_pre)


	rgb = np.zeros((4096, 4096, 3))
	rgb[:, :, 0] = ratio_a
	rgb[:, :, 1] = ratio_b
	rgb[:, :, 2] = ratio_c


	fig = plt.figure(figsize=(8, 8))
	ax = plt.subplot(projection=aia_171_map) 
	ax.imshow(rgb)
	aia_171_map.draw_grid(axes=ax,annotate=False)
##### Setting limit in both x and y axis
	xlims_world = [-350, 1225]*unit.arcsec
	ylims_world = [-1224, 150]*unit.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_171_map.coordinate_frame)
	pixel_coords = aia_171_map.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax.set_xlim(xlims_pixel)
	ax.set_ylim(ylims_pixel)
	ax.set_xlabel('X (arcsec)')
	ax.set_ylabel('Y (arcsec)')
	ax.set_title(aia_171_map.meta['date-obs'])
##### Saving each image
	print('saving image')
	plt.savefig(root+'aia_ratio_pyimages/image_{:03d}.png'.format(i))	
		 


