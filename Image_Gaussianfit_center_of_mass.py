##### This code is used to fit the radio images with Gaussian centroid and Center of mass (Choose a boundary of specific intensity and then mask it by 1 
##### and outside of that boundary is masked as zero).Then it calculates the shift of the sun in each image.
"""
==============================================================================
Create a Helioprojective Map from observations in the RA-DEC coordinate system
==============================================================================

How to create a `~sunpy.map.Map` in Helioprojective Coordinate Frame from radio observations
in GCRS (RA-DEC).

In this example a LOFAR FITS file (created with LOFAR's `Default Pre-Processing Pipeline (DPPP) and
WSClean Imager <https://support.astron.nl/LOFARImagingCookbook/dppp.html>`__) is read in,
the WCS header information is then used to make a new header with the information in Helioprojective,
and a `~sunpy.map.Map` is made.

The LOFAR example file has a WCS in celestial coordinates i.e. Right Ascension and
Declination (RA-DEC). For this example, we are assuming that the definition of LOFAR's
coordinate system for this observation is exactly the same as Astropy's ~astropy.coordinates.GCRS.
For many solar studies we may want to plot this data in some Sun-centered coordinate frame,
such as `~sunpy.coordinates.frames.Helioprojective`. In this example we read the data and
header information from the LOFAR FITS file and then create a new header with updated WCS
information to create a `~sunpy.map.Map` with a HPC coordinate frame. We will make use of the
`astropy.coordinates` and `sunpy.coordinates` submodules together with `~sunpy.map.make_fitswcs_header`
to create a new header and generate a `~sunpy.map.Map`.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
import pdb
import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import frames, sun
from scipy import optimize
from astropy.modeling.functional_models import Gaussian2D
from astropy.io import fits
from astropy.wcs import WCS
from pylab import *
from numpy import *    ### it is needed to define ravel
def fitgaussian(data):
     """Returns (height, x, y, width_x, width_y)
     the gaussian parameters of a 2D distribution found by a fit"""
     #params = moments(data)
     params =(200, 35, 30, 7, 10)
     errorfunction = lambda p: ravel(Gaussian2D(*p)(*indices(data.shape)) -data)
     p, success = optimize.leastsq(errorfunction, params) ## write verbose = 2 if you want to see the progress of the fitting
     return p
##############################################################################
###### We will first begin be reading in the header and data from the FITS file.
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576'
#root1 = '/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1-1095907872'

list1 = [root+'/QS_images/103-104/'+'1095907576_103-104_chan_53~55.fits',root + '/QS_images/093-094/'+'1095907576_093-094_chan_53~55.fits',root + '/QS_images/084-085/'+'1095907576_084-085_chan_53~55.fits',root + '/QS_images/076-077/'+'1095907576_076-077_chan_53~55.fits',root + '/QS_images/069-070/'+'1095907576_069-070_chan_53~55.fits',root + '/QS_images/062-063/'+'1095907576_062-063_chan_53~55.fits']


###### Plot image
fig = plt.figure(figsize=(8,6))
fig.subplots_adjust(hspace=0.4)
shift = np.zeros([6,2])
for i in range(len(list1)):
	hdu = fits.open(list1[i])[0]
	header = hdu.header
	data = hdu.data[0, 0, :, :]
	obstime = Time(header['date-obs'])
	frequency = header['crval3']*u.Hz
	lofar_loc = EarthLocation(lat=-26.7033*u.deg, lon=116.671*u.deg)
	lofar_gcrs = SkyCoord(lofar_loc.get_gcrs(obstime))
	reference_coord = SkyCoord(header['crval1']*u.Unit(header['cunit1']),
                           header['crval2']*u.Unit(header['cunit2']),
                           frame='gcrs',
                           obstime=obstime,
                           obsgeoloc=lofar_gcrs.cartesian,
                           obsgeovel=lofar_gcrs.velocity.to_cartesian(),
                           distance=lofar_gcrs.hcrs.distance)
	reference_coord_arcsec = reference_coord.transform_to(frames.Helioprojective(observer=lofar_gcrs))
	cdelt1 = (np.abs(header['cdelt1'])*u.deg).to(u.arcsec)
	cdelt2 = (np.abs(header['cdelt2'])*u.deg).to(u.arcsec)
	P1 = sun.P(obstime)
	new_header = sunpy.map.make_fitswcs_header(data, reference_coord_arcsec,
                                           reference_pixel=u.Quantity([header['crpix1']-1,
                                                                       header['crpix2']-1]*u.pixel),
                                           scale=u.Quantity([cdelt1, cdelt2]*u.arcsec/u.pix),
                                           rotation_angle=-P1,
                                           wavelength=frequency.to(u.MHz).round(2),
                                           observatory='LOFAR')
	lofar_map = sunpy.map.Map(data, new_header)
	lofar_map_rotate = lofar_map.rotate()
	bl = SkyCoord(-3000*u.arcsec, -3000*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
	tr = SkyCoord(3000*u.arcsec, 3000*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
	lofar_submap = lofar_map_rotate.submap(bl, top_right=tr)
	ax = fig.add_subplot(2,3,i+1,projection=lofar_submap)
	
	lofar_submap.plot(cmap='jet', 
        vmin=np.percentile(lofar_submap.data, 10.0),
        vmax=np.percentile(lofar_submap.data, 100.0))
	lofar_submap.draw_limb(color='black')
	lofar_submap.draw_grid(color='black')
	ax.set_xlabel(' ')
	ax.set_ylabel(' ')
	#lofar_submap.draw_contours(np.arange(50, 100, 5)*u.percent, colors='grey', linewidths=0.5)

	qs_boundary = (10.0/np.max(lofar_submap.data))*100.0
	lofar_submap.draw_contours([qs_boundary]*u.percent, colors='red', linewidths=0.8)

	####### Gaussian Fitting
	params = fitgaussian(lofar_submap.data)
	fit = Gaussian2D(*params)
	print(params)
	ax.contour(fit(*indices(lofar_submap.data.shape)), cmap=cm.copper)


	####### Plotting the centre of mass
	subdata = lofar_submap.data
	zeros = np.where(subdata < 10.0)
	ones = np.where(subdata > 10.0)
	subdata[zeros[0], zeros[1]] = 0.0
	subdata[ones[0], ones[1]] = 1.0
	import scipy.ndimage as ndi
	cy, cx = ndi.center_of_mass(subdata)
	ax.plot(cx,cy,'x', color='black')
	print(cx, cy)
	## Shift between the coordinates of the solar centre and the centroid. The pixel coordinates of the solar centre is (30.245, 30.4). 
	## Because the dimensions of the lofar_submap is 61 X 61, I start with 30, 30 for the coordinatrs of the solar centre.
	X_shift = (params[2]-30.245) *lofar_submap.meta['cdelt1']
	Y_shift = (params[1]-30.4) *lofar_submap.meta['cdelt1']
	shift[i][0] = X_shift /60
	shift[i][1] = Y_shift /60
	print(params[2]-cx) ##### the difference between the COM and Gaussian centroid position 
	print(params[1]-cy)
	print(X_shift/60)
	print(Y_shift/60)
	
#### Save the shift of the Sun's position for each image
#np.save(root +'/1095907576_103-104/shift_heliocoords.npy',shift)
plt.show()

