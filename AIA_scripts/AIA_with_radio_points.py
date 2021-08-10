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


def getarcsec(gauss_file, fint=0, tint=0):

	pos = np.load(gauss_file)

	xasec = pos[fint][tint][0]*u.arcsec
	yasec = pos[fint][tint][1]*u.arcsec

	return xasec.to('deg'), yasec.to('deg')

########Loading the AIA file

root = '/mnt/LOFAR-PSP/shilpi_MWA_project/'
aia_map = sunpy.map.Map(root+'AIA_data/aiaimages_171/aia_lev1_171a_2014_09_28t02_47_35_34z_image_lev1.fits')

############## Defining Axes

fig = plt.figure(figsize=(13, 7))
outer = gridspec.GridSpec(2, 3, wspace=0.5, hspace=0.2)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:,:1], wspace=0.1, hspace=0.3)

ax1 = plt.subplot(inner0[0], projection=aia_map)

sdoaia171 = plt.get_cmap('sdoaia171')
aia_map.plot(cmap = sdoaia171, norm=ImageNormalize(vmin=0,vmax=5000,stretch=LogStretch(10)))

##### Setting limit in both x and y axis
xlims_world = [-350, 1250]*unit.arcsec
ylims_world = [-2000, 200]*unit.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
pixel_coords = aia_map.world_to_pixel(world_coords)
# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value
ax1.set_xlim(xlims_pixel)
ax1.set_ylim(ylims_pixel)	
ax1.patch.set_facecolor('black')

from matplotlib import pylab
import matplotlib.colors as colors
cm = pylab.get_cmap('cool_r')
collist = cm(np.linspace(0, 255, 10).astype(int))



##############################################################################
##     Overplot the MWA points.
gauss_file = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')
xasec, yasec = getarcsec(gauss_file, fint=0, tint=173)
ax1.plot(xasec, yasec, 'o', color=collist[0], transform=ax1.get_transform('world'))



plt.show()

inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,1:], wspace=0.1, hspace=0.04)           ##### from 0th row and 1st column  

inner2 = gridspec.GridSpecFromSubplotSpec(6, 1,
        subplot_spec=outer[1:2,1:], wspace=0.1, hspace=0.05)
